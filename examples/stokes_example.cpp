/** @file stokes_example.cpp

    @brief Steady Stokes problem with adjoint approach for sensitivity
   analysis
*/

#include <gismo.h>

using namespace gismo;

// Global Typedefs
typedef gsExprAssembler<>::geometryMap geometryMap;
typedef gsExprAssembler<>::variable variable;
typedef gsExprAssembler<>::space space;
typedef gsExprAssembler<>::solution solution;

int main(int argc, char* argv[]) {
  ////////////////////
  // Global Options //
  ////////////////////

  // field IDs
  constexpr index_t PRESSURE_ID = 0;
  constexpr index_t VELOCITY_ID = 1;
  // field dimensions
  constexpr index_t PRESSURE_DIM = 1;
  // number of solution and test spaces
  constexpr index_t NUM_TRIAL = 2;
  constexpr index_t NUM_TEST = 2;

  // Setup values for timing
  double setup_time(0), assembly_time_ls(0), solving_time_ls(0),
      plotting_time(0);
  gsStopwatch timer;
  timer.restart();

  ////////////////////////////////
  // Parse Command Line Options //
  ////////////////////////////////

  // Title
  gsCmdLine cmd("Stokes Example");

  // Provide vtk data
  bool dont_plot = false;
  cmd.addSwitch("no-plot",
                "Suppress generation of a ParaView visualization file with " 
                "the solution", dont_plot);
  std::string output_file = "solutionStokes";
  cmd.addString("o", "outputfile", "Name of the ParaView output file.", output_file);
  bool export_xml = false;
  cmd.addSwitch("export-xml", "Export solution into g+smo xml format.",
                export_xml);
  index_t sample_rate{4};
  cmd.addInt("q", "sample-rate", "Sample rate of splines for export",
             sample_rate);

  // Material constants
  real_t viscosity{6000};
  cmd.addReal("v", "visc", "Viscosity", viscosity);

  // Mesh options
  index_t degreeElevate = 0;
  cmd.addInt("e", "degreeElevate", "Number of uniform degree elevations",
             degreeElevate);
  index_t numRefine = 0;
  cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement loops",
             numRefine);

  std::string fn("../../filedata/pde/microstructure.xml");
  cmd.addString("f", "file", "Input XML file", fn);

  // A few more mesh options
  index_t mp_id{0}, vel_bc_id{2}, p_bc_id{3}, ass_opt_id{4};
  cmd.addInt("m", "multipach_id", "ID of the multipatch mesh in mesh file",
             mp_id);
  cmd.addInt("b", "boundary_id",
             "ID of the boundary condition function in mesh file", vel_bc_id);
  cmd.addInt("a", "assembly_options_id",
             "ID of the assembler options in mesh file", ass_opt_id);
#ifdef _OPENMP
  int numThreadsRequested{1};
  cmd.addInt("p", "n_threads", "Number of threads used", numThreadsRequested);
#endif

  // Parse command line options
  try {
    cmd.getValues(argc, argv);
  } catch (int rv) {
    return rv;
  }

  // Import mesh and load relevant information
  gsFileData<> fd(fn);
  gsInfo << "Loaded file " << fd.lastPath() << std::endl;

  // retrieve multi-patch data
  gsMultiPatch<> domain_patches;
  fd.getId(mp_id, domain_patches);

  // retrieve velocity boundary conditions from file
  gsBoundaryConditions<> velocity_bcs;
  fd.getId(vel_bc_id, velocity_bcs);
  velocity_bcs.setGeoMap(domain_patches);
  gsInfo << "Velocity boundary conditions:\n" << velocity_bcs << std::endl;

  bool has_p_bc = false;
  gsBoundaryConditions<> pressure_bcs;
  if (fd.hasId(p_bc_id)) {
    has_p_bc = true;
    fd.getId(p_bc_id, pressure_bcs);
    pressure_bcs.setGeoMap(domain_patches);
    gsInfo << "Presure boundary conditions:\n" << pressure_bcs << '\n';
  }

  // retrieve assembly options
  gsOptionList Aopt;
  fd.getId(ass_opt_id, Aopt);

  const index_t geomDim = domain_patches.geoDim();
  gsInfo << "Geometric dimension " << geomDim << std::endl;

  //! Create function bases for pressure and velocity
  gsMultiBasis<> function_basis_pressure(
      domain_patches,
      true);  // true: poly-splines (not NURBS)
  gsMultiBasis<> function_basis_velocity(
      domain_patches,
      true);  // true: poly-splines (not NURBS)

  // Elevate the degree as desired by the user and increase the degree one 
  // additional time for the velocity to obtain Taylor-Hood elements
  function_basis_pressure.setDegree( 
      function_basis_pressure.maxCwiseDegree() + degreeElevate
  );
  function_basis_velocity.setDegree( 
      function_basis_velocity.maxCwiseDegree() + degreeElevate + 1
  );

  // h-refine each basis (for performing the analysis)
  for (int r=0; r<numRefine; ++r) {
    function_basis_pressure.uniformRefine();
    function_basis_velocity.uniformRefine();
  }

  // Output user information
  gsInfo << "Summary Velocity:" << std::endl
         << "Patches: " << domain_patches.nPatches()
         << "\nMin-degree: " << function_basis_velocity.minCwiseDegree() 
         << "\nMax-degree: " << function_basis_velocity.maxCwiseDegree()
         << std::endl << std::endl;
  gsInfo << "Summary Pressure:" << std::endl
         << "Patches: " << domain_patches.nPatches()
         << "\nMin-degree: " << function_basis_pressure.minCwiseDegree()
         << "\nMax-degree: " << function_basis_pressure.maxCwiseDegree()
         << std::endl << std::endl;
#ifdef _OPENMP
  index_t maxOmpThreads = omp_get_max_threads();
  index_t numThreadsUsed = std::min(maxOmpThreads, numThreadsRequested);
  gsInfo << "Available threads: " << maxOmpThreads << std::endl;
  gsInfo << "Number of threads: " << numThreadsUsed << std::endl << std::endl;
  omp_set_num_threads(numThreadsUsed);
#endif

  ///////////////////
  // Problem Setup //
  ///////////////////

  // Construct expression assembler
  // (takes number of test and solution function spaces as arguments)
  gsExprAssembler<> expr_assembler(NUM_TEST, NUM_TRIAL);
  expr_assembler.setOptions(Aopt);
  gsInfo << "Active options:\n" << expr_assembler.options() << std::endl;

  // Elements used for numerical integration
  expr_assembler.setIntegrationElements(function_basis_velocity);

  // Set the geometry map
  geometryMap geoMap = expr_assembler.getMap(domain_patches);

  // Set the discretization spaces
  space pressure_trial_space = expr_assembler.getSpace(
    function_basis_pressure, PRESSURE_DIM, PRESSURE_ID
  );
  gsInfo << "Solution space for pressure (id=" << pressure_trial_space.id() 
         << ") has " << pressure_trial_space.rows() << " rows and " 
         << pressure_trial_space.cols() << " columns." << std::endl;
  
  space velocity_trial_space = expr_assembler.getSpace(
    function_basis_velocity, geomDim, VELOCITY_ID
  );
  gsInfo << "Solution space for velocity (id=" << velocity_trial_space.id() 
         << ") has " << velocity_trial_space.rows() << " rows and " 
         << velocity_trial_space.cols() << " columns." << std::endl;

  // Solution vector and solution variable
  gsMatrix<> full_solution;
  solution pressure_field = 
      expr_assembler.getSolution(pressure_trial_space, full_solution);
  solution velocity_field =
      expr_assembler.getSolution(velocity_trial_space, full_solution);

  // Intitalize multi-patch interfaces for pressure field
  if (has_p_bc) {
    pressure_trial_space.setup(pressure_bcs, dirichlet::l2Projection, 0);
  } else {
    pressure_trial_space.setup(0);
  }
  // Initialize interfaces and Dirichlet bcs for velocity field
  velocity_trial_space.setup(velocity_bcs, dirichlet::l2Projection, 0);

  // Initialize the system
  expr_assembler.initSystem();
  setup_time += timer.stop();

  gsInfo << "Number of degrees of freedom : " << expr_assembler.numDofs()
         << std::endl;
  gsInfo << "Number of blocks in the system matrix : "
         << expr_assembler.numBlocks() << std::endl;

  //////////////
  // Assembly //
  //////////////

  gsInfo << "Starting assembly of linear system ..." << std::flush;
  timer.restart();

  // Compute the system matrix and right-hand side
  // Variational formulation
  // ∫ q tr(∇v) + μ ∇v:∇w + μ (∇v)ᵀ:∇w - p tr(∇w) dV = 0
  auto phys_jacobian = ijac(velocity_trial_space, geoMap);                                            // ∇v
  auto bilin_conti = pressure_trial_space * idiv(velocity_trial_space, geoMap).tr() * meas(geoMap);
  auto bilin_press = -idiv(velocity_trial_space, geoMap) * pressure_trial_space.tr() * meas(geoMap);
  auto bilin_mu_1 = viscosity * (phys_jacobian.cwisetr() % phys_jacobian.tr()) *
                    meas(geoMap);
  auto bilin_mu_2 =
      viscosity * (phys_jacobian % phys_jacobian.tr()) * meas(geoMap);

  expr_assembler.assemble(bilin_conti, bilin_press, bilin_mu_1, bilin_mu_2);

  assembly_time_ls += timer.stop();
  gsInfo << "\t\tFinished" << std::endl;

  ///////////////////
  // Linear Solver //
  ///////////////////

  gsInfo << "Solving the linear system of equations ..." << std::flush;
  timer.restart();

  const auto& system_matrix = expr_assembler.matrix();
  const auto& rhs_vector = expr_assembler.rhs();

  // Initialize linear solver
  gsSparseSolver<>::BiCGSTABILUT solver;
//   gsSparseSolver<>::LeastSquaresCG solver;
  solver.compute(system_matrix);
  full_solution = solver.solve(rhs_vector);

  solving_time_ls += timer.stop();
  gsInfo << "\tFinished" << std::endl;

  ////////////////////
  // Postprocessing //
  ////////////////////

  gsExprEvaluator<> expression_evaluator(expr_assembler);

  //////////////////////////////
  // Export and Visualization //
  //////////////////////////////

  // Generate Paraview File
  if (!dont_plot) {
    gsInfo << "\nStarting the paraview export ..." << std::flush;
    timer.restart();

    gsParaviewCollection collection("ParaviewOutput/" + output_file,
                                    &expression_evaluator);
    collection.options().setSwitch("plotElements", false);
    collection.options().setInt("plotElements.resolution", sample_rate);
    collection.newTimeStep(&domain_patches);
    collection.addField(pressure_field, "pressure");
    collection.addField(velocity_field, "velocity");
    collection.saveTimeStep();
    collection.save();

    plotting_time += timer.stop();
    gsInfo << "\tFinished" << std::endl;
  }

  // Export solution file as xml
  if (export_xml) {
    gsInfo << "Starting the xml export ..." << std::flush;

    // Export pressure
    gsMatrix<> full_solution_pressure;
    gsFileData<> output_pressure;
    pressure_field.extractFull(
        full_solution_pressure);  // patch-wise solution with BCs
    output_pressure << full_solution_pressure;
    output_pressure.save("pressure_field.xml");

    // Export velocity
    gsMatrix<> full_solution_velocity;
    gsFileData<> output;
    velocity_field.extractFull(
        full_solution_velocity);  // patch-wise solution with BCs
    output << full_solution_velocity;
    output.save("velocity_field.xml");

    gsInfo << "\t\tFinished" << std::endl;
  }

  // User output infor timings
  gsInfo << "\n\nTotal time: "
         << setup_time + assembly_time_ls + solving_time_ls + plotting_time
         << std::endl;
  gsInfo << "                       Setup: " << setup_time << std::endl;
  gsInfo << "      Assembly Linear System: " << assembly_time_ls << std::endl;
  gsInfo << "       Solving Linear System: " << solving_time_ls << std::endl;
  gsInfo << "                    Plotting: " << plotting_time << std::endl
         << std::flush;

  return EXIT_SUCCESS;

}  // end main