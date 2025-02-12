/* @file stokesConvergenceStudy.cpp

  @brief Stationary Stokes with options of h- and p-refinement.
  Should serve as test if running the same simulation multiple times with uneven degree of velocity
  sometimes results in different solutions. The main culprit is the different rhs'es. The 
  system matrix seems to be the same for everything. It might be due to the imposition of the
  boundary conditions.
*/

#include <gismo.h>

using namespace gismo;

// Global Typedefs
typedef gsExprAssembler<>::geometryMap geometryMap;
typedef gsExprAssembler<>::variable variable;
typedef gsExprAssembler<>::space space;
typedef gsExprAssembler<>::solution solution;

template<typename T>
void printDim(T expression, std::string name) {
  gsInfo << "Dimension of " << name << ": (" << expression.rows() << "," << expression.cols() << ")\n";
}

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
  gsCmdLine cmd("Stokes Example - h-/p-convergence analysis");

  // Provide vtk data
  bool dont_plot = false;
  cmd.addSwitch("no-plot",
                "Suppress generation of a ParaView visualization file with " 
                "the solution", dont_plot);
  std::string output_file = "solutionStokesConvergence_Vector";
  cmd.addString("o", "outputfile", "Name of the ParaView output file.", output_file);
  index_t sample_rate{4};
  bool output_lhs_rhs = false;
  cmd.addSwitch("output-lhs-rhs",
                "Output system matrix and rhs vector in command line",
                output_lhs_rhs);
  bool zero_pressure = false;
  cmd.addSwitch("zero-pressure",
                "Set integral of pressure over domain to zero and disable pressure BCs",
                zero_pressure);
  bool fix_pressure_at_se_corner = false;
  cmd.addSwitch("pressure-se-corner",
                "Fix pressure at southeast corner of patch 0 to 0",
                fix_pressure_at_se_corner
  );
  // cmd.addInt("q", "sample-rate", "Sample rate of splines for export",
  //            sample_rate);

  // Material constants
  real_t viscosity{1};
  cmd.addReal("v", "visc", "Viscosity", viscosity);

  // Mesh options
  index_t hRef{0}, pRef{0};
  cmd.addInt("r", "href", "Number of uniform h-refinements", hRef);
  cmd.addInt("p", "pref", "Number of p-refinements", pRef);


  std::string fn("../../filedata/stokesConvergenceStudy/square.xml");
  cmd.addString("f", "file", "Input XML file", fn);

  // A few more mesh options
  index_t mp_id{0}, vel_bc_id{1}, ass_opt_id{10}, p_bc_id{2}, body_force_id{100},
        vel_analytical_id{101}, p_analytical_id{102};
  cmd.addInt("m", "multipach_id", "ID of the multipatch mesh in mesh file",
             mp_id);
  cmd.addInt("b", "boundary_id",
             "ID of the velocity boundary condition function in mesh file", vel_bc_id);
  cmd.addInt("a", "assembly_options_id",
             "ID of the assembler options in mesh file", ass_opt_id);
  // cmd.addInt("p", "p_bc_id",
  //            "ID of the assembler options in mesh file", p_bc_id);
// #ifdef _OPENMP
//   int numThreadsRequested{1};
//   cmd.addInt("p", "n_threads", "Number of threads used", numThreadsRequested);
// #endif

  // Parse command line options
  try {
    cmd.getValues(argc, argv);
  } catch (int rv) {
    return rv;
  }

  if (zero_pressure && fix_pressure_at_se_corner) {
    GISMO_ERROR("Cannot have both zero pressure integral and fixing pressure at corner. Choose either one of those options or none.");
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

  // retrieve pressure boundary conditions from file
  gsBoundaryConditions<> pressure_bcs;

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
      function_basis_pressure.maxCwiseDegree() + pRef
  );
  function_basis_velocity.setDegree( 
      function_basis_velocity.maxCwiseDegree() + pRef + 1
  );

  // h-refine each basis (for performing the analysis)
  for (int r=0; r<hRef; ++r) {
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
  gsExprEvaluator<> expression_evaluator(expr_assembler);

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

  // gsSparseSolver<>::BiCGSTABILUT solver;
//   gsSparseSolver<>::LeastSquaresCG solver;
  gsSparseSolver<>::LU solver;
  gsSparseSolver<>::BiCGSTABILUT solver_iterative;


  gsFunctionExpr<> body_force_func;
  fd.getId(body_force_id, body_force_func);
  auto body_force = expr_assembler.getCoeff(body_force_func, geoMap);

  // gismo::dirichlet::values const &l2Projection = gismo::dirichlet::l2Projection;

  // Intitalize multi-patch interfaces for pressure field
  if (!zero_pressure) {
    if (fix_pressure_at_se_corner) {
      gsInfo << "Ignoring pressure BCs from file. Instead setting the pressure at the southeast corner of patch 0 to 0.\n";
      gismo::gsConstantFunction<> const zero(0.0, domain_patches.targetDim());
      pressure_bcs.addCondition(0, boundary::southeast, condition_type::dirichlet, zero, PRESSURE_ID);
      pressure_bcs.setGeoMap(domain_patches);
      pressure_trial_space.setup(pressure_bcs, dirichlet::l2Projection, 0);
    } else {
      // Else apply presure BCs from xml file
      fd.getId(p_bc_id, pressure_bcs);
      pressure_bcs.setGeoMap(domain_patches);
      pressure_trial_space.setup(pressure_bcs, Aopt.getInt("DirichletValues"), 0);
    }
  }
  gsInfo << "Pressure boundary conditions:\n" << pressure_bcs << std::endl;
  // Initialize interfaces and Dirichlet bcs for velocity field
  velocity_trial_space.setup(velocity_bcs, Aopt.getInt("DirichletValues"), 0);
  // pressure_trial_space.setup(pressure_bcs, l2Projection, 0);
  // velocity_trial_space.setup(velocity_bcs, l2Projection, 0);

  // Initialize the system
  expr_assembler.initSystem();
  setup_time += timer.stop();

  gsInfo << expr_assembler.numDofs();
  // gsInfo << "Number of degrees of freedom : " << expr_assembler.numDofs()
  //       << std::endl;
  // gsInfo << "Number of blocks in the system matrix : "
  //       << expr_assembler.numBlocks() << std::endl;

  //////////////
  // Assembly //
  //////////////

  // Starting assembly of linear system
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
  // auto ext_force = (body_force.cwisetr() * velocity_trial_space.tr()) * meas(geoMap);
  // auto ext_force = (velocity_trial_space % body_force.tr()) * meas(geoMap);
  // auto ext_force1 = body_force.cwisetr()[0] * velocity_trial_space[0].tr() * meas(geoMap);
  // auto ext_force2 = body_force.cwisetr()[1] * velocity_trial_space[1].tr() * meas(geoMap);
  auto ext_force1 = velocity_trial_space[0] * body_force.cwisetr()[0].tr() * meas(geoMap);
  auto ext_force2 = velocity_trial_space[1] * body_force.cwisetr()[1].tr() * meas(geoMap);

  auto zero_integral_pressure = pressure_trial_space * pressure_trial_space.tr() * meas(geoMap);

  if (!zero_pressure) {
    expr_assembler.assemble(bilin_conti, bilin_press, bilin_mu_1, bilin_mu_2, ext_force1, ext_force2);
  } else {
    expr_assembler.assemble(bilin_conti, bilin_press, bilin_mu_1, bilin_mu_2, ext_force1, ext_force2, zero_integral_pressure);
  }

  assembly_time_ls += timer.stop();
  gsInfo << "." << std::flush;

  ///////////////////
  // Linear Solver //
  ///////////////////

  // gsInfo << "Solving the linear system of equations ..." << std::flush;
  timer.restart();

  const auto& system_matrix = expr_assembler.matrix();
  const auto& rhs_vector = expr_assembler.rhs();

  // Initialize linear solver
  solver.compute(system_matrix);
  full_solution = solver.solve(rhs_vector);

  // solver_iterative.setTolerance(1e-100);
  // solver_iterative.setMaxIterations(1000);
  // solver_iterative.compute(system_matrix);
  // // full_solution = solver_iterative.solveWithGuess(rhs_vector, full_solution);
  // full_solution = solver_iterative.solve(rhs_vector);
  // gsInfo << "----------------------------------------------\n";
  // gsInfo << "Iterative error: " << solver_iterative.error() << '\n';
  // gsInfo << "Iterations: " << solver_iterative.iterations() << '\n';
  // // gsInfo << "Tolerance: " << solver_iterative.tolerance() << '\n';

  // gsMatrix<> iterative_error_history;
  // solver_iterative.initIteration(rhs_vector, full_solution);
  // solver_iterative.solveDetailed(rhs_vector, full_solution, iterative_error_history);

  solving_time_ls += timer.stop();
  gsInfo << "." << std::flush;

  // Error evaluation
  gsFunctionExpr<> vel_analytical_func, p_analytical_func;
  fd.getId(vel_analytical_id, vel_analytical_func);
  fd.getId(p_analytical_id, p_analytical_func);
  auto vel_analytical = expression_evaluator.getVariable(vel_analytical_func, geoMap);
  auto p_analytical = expression_evaluator.getVariable(p_analytical_func, geoMap);

  real_t velocity_error, pressure_error;
  velocity_error = math::sqrt(expression_evaluator.integral(
    (vel_analytical - velocity_field).sqNorm() * meas(geoMap)
  ));
  pressure_error = math::sqrt(expression_evaluator.integral(
    (p_analytical - pressure_field).sqNorm() * meas(geoMap)
  ));


  gsInfo << "Errors: " << "\n"
    << "Velocity: " << velocity_error << '\n'
    << "Pressure: " << pressure_error << '\n';

  

  ////////////////////
  // Postprocessing //
  ////////////////////

  if (output_lhs_rhs) {
    gsInfo << "lhs---\n" << expr_assembler.matrix() << "\n---lhs";
    gsInfo << "rhs---\n" << expr_assembler.rhs()    << "\n---rhs";
    gsInfo << "sol---\n" << full_solution           << "\n---sol";
  }

  //////////////////////////////
  // Export and Visualization //
  //////////////////////////////

  // Generate Paraview File
  if (!dont_plot) {
    gsInfo << "\nStarting the paraview export ..." << std::flush;
    timer.restart();

    gsParaviewCollection collection("ParaviewOutput/" + output_file,
                                    &expression_evaluator);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", sample_rate);
    collection.options().setInt("precision", 12);
    collection.newTimeStep(&domain_patches);
    collection.addField(pressure_field, "pressure");
    collection.addField(velocity_field, "velocity");
    collection.addField(vel_analytical, "velocity (analytical)");
    collection.addField(p_analytical, "pressure (analytical)");
    collection.addField(body_force, "Body force");
    collection.addField(pressure_field - p_analytical, "Pressure error");
    collection.addField(velocity_field - vel_analytical, "Velocity error");
    collection.saveTimeStep();
    collection.save();

    plotting_time += timer.stop();
    gsInfo << "\tFinished" << std::endl;
  }

  // User output infor timings
  gsInfo << "\n\nTotal time: "
         << setup_time + assembly_time_ls + solving_time_ls + plotting_time
         << std::endl;

  return EXIT_SUCCESS;

}  // end main