/** @file stokes_recreation.cpp

    @brief Recreation and computation of recreated velocity and pressure field for 
    stationary Stokes equation using microstructured geometry
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

  ////////////////////////////////
  // Parse Command Line Options //
  ////////////////////////////////

  // Title
  gsCmdLine cmd("Stokes Recreation Example");

  // Provide vtk data
  bool dont_plot = false;
  cmd.addSwitch("no-plot",
                "Suppress generation of a ParaView visualization file with " 
                "the solution", dont_plot);
  std::string output_file = "recreationStokes";
  cmd.addString("o", "outputfile", "Name of the ParaView output file.", output_file);
  index_t sample_rate{4};
  cmd.addInt("s", "sample-rate", "Sample rate of splines for export",
             sample_rate);

  // Mesh options
  index_t degreeElevate = 0;
  cmd.addInt("e", "degreeElevate", "Number of uniform degree elevations",
             degreeElevate);
  index_t numRefine = 0;
  cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement loops",
             numRefine);

  // Input files
  std::string vf("velocity_field.xml"), pf("pressure_field.xml"),
              vfr("velocity_field_recreation.xml"),
              pfr("pressure_field_recreation.xml"),
              fn("microstructure.xml");
  cmd.addString("v", "velocity-file", "Input file for velocity field", vf);
  cmd.addString("p", "pressure-file", "Input file for pressure field", pf);
  cmd.addString("w", "velocity-recreation", "Input file of recreated velocity field", vfr);
  cmd.addString("q", "pressure-recreation", "Input file of recreated pressure field", pfr);
  cmd.addString("f", "microstructure-file", "Input file for geometry", fn);

  // A few more mesh options
  index_t mp_id{0};
  cmd.addInt("m", "multipach_id", "ID of the multipatch mesh in mesh file",
             mp_id);

  // Parse command line options
  try {
    cmd.getValues(argc, argv);
  } catch (int rv) {
    return rv;
  }

  // Import mesh and load relevant information
  gsFileData<> fd(fn);
  gsInfo << "Loaded geometry file: " << fd.lastPath() << std::endl;

  // retrieve multi-patch data
  gsMultiPatch<> domain_patches;
  fd.getId(mp_id, domain_patches);

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

  ///////////////////
  // Problem Setup //
  ///////////////////

  // Construct expression assembler
  // (takes number of test and solution function spaces as arguments)
  gsExprAssembler<> expr_assembler(NUM_TEST, NUM_TRIAL);
  // gsInfo << "Active options:\n" << expr_assembler.options() << std::endl;

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

  // -------------------------------------------------------------------
  // Read the solution vectors
  gsFileData<> vfd(vf), pfd(pf), vrfd(vfr), prfd(pfr);
  gsInfo << "Read the solution vector xml files.\n";

  gsMatrix<> solution_pressure, solution_velocity, solution_pressure_recreated,
            solution_velocity_recreated;
  vfd.getId(0, solution_velocity);
  pfd.getId(0, solution_pressure);
  vrfd.getId(0, solution_velocity_recreated);
  prfd.getId(0, solution_pressure_recreated);
  gsInfo << "Sizes of solution vectors:\n"
      << "\tVelocity: " << solution_velocity.dim() << '\n'
      << "\tPressure: " << solution_pressure.dim() << '\n';

  gsInfo << "Sizes of bases:\n"
      << "\tVelocity: " << function_basis_velocity.size() << '\n'
      << "\tPressure: " << function_basis_pressure.size() << '\n';


  // Solution vector and solution variable
  solution pressure_field = 
      expr_assembler.getSolution(pressure_trial_space, solution_pressure);
  solution velocity_field =
      expr_assembler.getSolution(velocity_trial_space, solution_velocity);
  solution pressure_field_recreated =
      expr_assembler.getSolution(pressure_trial_space, solution_pressure_recreated);
  solution velocity_recreated = 
      expr_assembler.getSolution(velocity_trial_space, solution_velocity_recreated);

  // Initialize the system
  expr_assembler.initSystem();

  gsInfo << "Number of degrees of freedom : " << expr_assembler.numDofs()
         << std::endl;

  ////////////////////
  // Postprocessing //
  ////////////////////

  gsExprEvaluator<> expression_evaluator(expr_assembler);

  // gsInfo << "\n-------------------------\nErrors:\n";
  gsInfo << "   ___  ____   ____   ___   ____    _____    \n" <<
"  /  _]|    \\ |    \\ /   \\ |    \\  / ___/ __ \n" <<
" /  [_ |  D  )|  D  )     ||  D  )(   \\_ |  |\n" <<
"|    _]|    / |    /|  O  ||    /  \\__  ||__|\n" <<
"|   [_ |    \\ |    \\|     ||    \\  /  \\ | __ \n" <<
"|     ||  .  \\|  .  \\     ||  .  \\ \\    ||  |\n" <<
"|_____||__|\\_||__|\\_|\\___/ |__|\\_|  \\___||__|\n";

  gsMatrix<> pressure_difference = solution_pressure - solution_pressure_recreated;
  solution error_pressure =
      expr_assembler.getSolution(pressure_trial_space, pressure_difference);

  double l1error_pressure = expression_evaluator.integral(
    abs(pressure_field.val() - pressure_field_recreated.val()) * meas(geoMap)
  );
  double l1error_rel_pressure = expression_evaluator.integral(
    (abs((pressure_field.val() - pressure_field_recreated.val()) / pressure_field.val()))
    * meas(geoMap)
  );
  double l2error_pressure = math::sqrt(expression_evaluator.integral(
    (pressure_field - pressure_field_recreated).sqNorm() * meas(geoMap)
  ));
  double l2error_rel_pressure = math::sqrt(expression_evaluator.integral(
    ((pressure_field - pressure_field_recreated) / pressure_field.val()).sqNorm()
    * meas(geoMap)
  ));
  double linferror_pressure = expression_evaluator.max(
    abs(pressure_field.val() - pressure_field_recreated.val())
  );
  double linferror_rel_pressure = expression_evaluator.max(
    abs(pressure_field.val() - pressure_field_recreated.val()) / abs(pressure_field.val())
  );
  double h1error_pressure = l2error_pressure + math::sqrt(expression_evaluator.integral(
    (igrad(pressure_field) - igrad(pressure_field_recreated)).sqNorm() * meas(geoMap)
  ));

  gsInfo << std::setprecision(16);
  gsInfo  << "Pressure:"
          << "\n\tL1        : " << l1error_pressure
          << "\n\tL2        : " << l2error_pressure
          << "\n\tLinf      : " << linferror_pressure
          << "\n\tH1        : " << h1error_pressure
          << "\n\tL1  (rel.): " << l1error_rel_pressure
          << "\n\tL2  (rel.): " << l2error_rel_pressure
          << "\n\tLinf(rel.): " << linferror_rel_pressure
          << '\n';

  //////////////////////////////
  // Export and Visualization //
  //////////////////////////////

  // Generate Paraview File
  if (!dont_plot) {
    gsInfo << "\nStarting the paraview export ..." << std::flush;

    gsParaviewCollection collection("ParaviewOutput/" + output_file,
                                    &expression_evaluator);
    collection.options().setSwitch("plotElements", false);
    collection.options().setInt("plotElements.resolution", sample_rate);
    collection.newTimeStep(&domain_patches);
    collection.addField(pressure_field, "pressure");
    collection.addField(pressure_field_recreated, "pressure (recreated)");
    collection.addField(error_pressure.norm(), "pressure (Error (abs.))");
    collection.addField(
      error_pressure.norm() / pressure_field.val(),
      "pressure (Error (rel.))"
    );
    // collection.addField(velocity_field, "velocity");
    collection.saveTimeStep();
    collection.save();
    gsInfo << "\tFinished" << std::endl;
  }

  return EXIT_SUCCESS;

}  // end main