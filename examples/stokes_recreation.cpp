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
  index_t mp_id{0}, vel_bc_id{2}, p_bc_id{3};
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

  GISMO_ASSERT(solution_velocity.size() == solution_velocity_recreated.size(), 
    "Size of solution vectors for velocity do not have the same dimension!");
  GISMO_ASSERT(solution_pressure.size() == solution_pressure_recreated.size(), 
    "Size of solution vectors for pressure do not have the same dimension!");

  index_t psize = solution_pressure.size();
  index_t vsize = solution_velocity.size();
  gsMatrix<> full_solution(psize+vsize, 1), full_solution_recreated(psize+vsize, 1);
  for (index_t i = 0; i < psize; i++) {
    full_solution[i] = solution_pressure[i];
    full_solution_recreated[i] = solution_pressure_recreated[i];
  }
  for (index_t i = 0; i < vsize; i++) {
    full_solution[psize+i] = solution_velocity[i];
    full_solution_recreated[psize+i] = solution_velocity_recreated[i];
  }

  // Solution vector and solution variable
  solution pressure_field = 
      expr_assembler.getSolution(pressure_trial_space, full_solution);
  solution velocity_field =
      expr_assembler.getSolution(velocity_trial_space, full_solution);
  solution pressure_field_recreated =
      expr_assembler.getSolution(pressure_trial_space, full_solution_recreated);
  solution velocity_field_recreated = 
      expr_assembler.getSolution(velocity_trial_space, full_solution_recreated);

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

  gsInfo << std::setprecision(24);
  gsInfo  << "Pressure:"
          << "\n\tL1        : " << l1error_pressure
          << "\n\tL2        : " << l2error_pressure
          << "\n\tLinf      : " << linferror_pressure
          << "\n\tH1        : " << h1error_pressure
          << "\n\tL1  (rel.): " << l1error_rel_pressure
          << "\n\tL2  (rel.): " << l2error_rel_pressure
          << "\n\tLinf(rel.): " << linferror_rel_pressure
          << '\n';

  // Workaround for velocity: get own solution field and only set respective entries
  // in the solution vector
  gsMatrix<> full_velx_solution(full_solution.size(), 1), full_vely_solution(full_solution.size(), 1),
            full_velx_rec_solution(full_solution.size(), 1), full_vely_rec_solution(full_solution.size(), 1);
  for (index_t i = 0; i < vsize; i++) {
    index_t xIndex = i+psize;
    index_t yIndex = i + psize+vsize/2;
    full_velx_solution[xIndex] = full_solution[xIndex];
    full_velx_rec_solution[xIndex] = full_solution_recreated[xIndex];
    full_vely_solution[yIndex] = full_solution[yIndex];
    full_vely_rec_solution[yIndex] = full_solution_recreated[yIndex];
  }

  solution velx_field = expr_assembler.getSolution(velocity_trial_space, full_velx_solution);
  solution velx_rec_field = expr_assembler.getSolution(velocity_trial_space, full_velx_rec_solution);
  solution vely_field = expr_assembler.getSolution(velocity_trial_space, full_vely_solution);
  solution vely_rec_field = expr_assembler.getSolution(velocity_trial_space, full_vely_rec_solution);
  
  gsMatrix<> velx_difference = full_velx_solution - full_velx_rec_solution;
  gsMatrix<> vely_difference = full_vely_solution - full_vely_rec_solution;
  solution error_velx = expr_assembler.getSolution(velocity_trial_space, velx_difference);
  solution error_vely = expr_assembler.getSolution(velocity_trial_space, vely_difference);

  gsInfo << "After the error getSolutions\n";

  double l1error_velx = expression_evaluator.integral(
    abs(velx_field.val() - velx_rec_field.val()) * meas(geoMap)
  );
  gsInfo << "After first error calculation\n";
  double l1error_rel_velx = expression_evaluator.integral(
    (abs(velx_field.val() - velx_rec_field.val()) / velx_field.val()) * meas(geoMap)
  );
  double l2error_velx = math::sqrt(expression_evaluator.integral(
    (velx_field - velx_rec_field).sqNorm() * meas(geoMap)
  ));
  double l2error_rel_velx = math::sqrt(expression_evaluator.integral(
    ((velx_field - velx_rec_field) / velx_field.val()).sqNorm() * meas(geoMap)
  ));
  double linferror_velx = expression_evaluator.max(
    abs(velx_field.val() - velx_rec_field.val())
  );
  double linferror_rel_velx = expression_evaluator.max(
    abs(velx_field.val() - velx_rec_field.val()) / abs(velx_field.val())
  );
  double h1error_velx = l2error_velx + math::sqrt(expression_evaluator.integral(
    (igrad(velx_field) - igrad(velx_rec_field)).sqNorm() * meas(geoMap)
  ));

  gsInfo  << "Velocity (x):"
          << "\n\tL1        : " << l1error_velx
          << "\n\tL2        : " << l2error_velx
          << "\n\tLinf      : " << linferror_velx
          << "\n\tH1        : " << h1error_velx
          << "\n\tL1  (rel.): " << l1error_rel_velx
          << "\n\tL2  (rel.): " << l2error_rel_velx
          << "\n\tLinf(rel.): " << linferror_rel_velx
          << '\n';

  double l1error_vely = expression_evaluator.integral(
    abs(vely_field.val() - vely_rec_field.val()) * meas(geoMap)
  );
  // double l1error_vely = expression_evaluator.integral(
  //   (vely_field.val() - vely_rec_field.val()) * (vely_field.val() - vely_rec_field.val()) * meas(geoMap)
  // );
  double l1error_rel_vely = expression_evaluator.integral(
    (abs(vely_field.val() - vely_rec_field.val()) / vely_field.val()) * meas(geoMap)
  );
  double l2error_vely = math::sqrt(expression_evaluator.integral(
    (vely_field - vely_rec_field).sqNorm() * meas(geoMap)
  ));
  double l2error_rel_vely = math::sqrt(expression_evaluator.integral(
    ((vely_field - vely_rec_field) / vely_field.val()).sqNorm() * meas(geoMap)
  ));
  double linferror_vely = expression_evaluator.max(
    abs(vely_field.val() - vely_rec_field.val())
  );
  double linferror_rel_vely = expression_evaluator.max(
    abs(vely_field.val() - vely_rec_field.val()) / abs(vely_field.val())
  );
  double h1error_vely = l2error_vely + math::sqrt(expression_evaluator.integral(
    (igrad(vely_field) - igrad(vely_rec_field)).sqNorm() * meas(geoMap)
  ));

  gsInfo  << "Velocity (y):"
          << "\n\tL1        : " << l1error_vely
          // << "\n\tyolo      : " << expression_evaluator.integral(velx_rec_field.norm() * meas(geoMap))
          // << "\n\tpressure  : " << expression_evaluator.integral(pressure_field_recreated * meas(geoMap))
          << "\n\tL2        : " << l2error_vely
          << "\n\tLinf      : " << linferror_vely
          << "\n\tH1        : " << h1error_vely
          << "\n\tL1  (rel.): " << l1error_rel_vely
          << "\n\tL2  (rel.): " << l2error_rel_vely
          << "\n\tLinf(rel.): " << linferror_rel_vely
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
    collection.addField(velocity_field, "velocity");
    collection.addField(velocity_field_recreated, "velocity (recreated)");
    collection.addField(error_velx.norm(), "velocity_x (Error (abs.))");
    collection.addField(error_velx.norm() / velx_field.val(), "velocity_x (Error (rel.))");
    collection.addField(error_vely.norm(), "velocity_y (Error (abs.))");
    collection.addField(error_vely.norm() / vely_field.val(), "velocity_y (Error (rel.))");
    collection.saveTimeStep();
    collection.save();
    gsInfo << "\tFinished" << std::endl;
  }

  return EXIT_SUCCESS;

}  // end main