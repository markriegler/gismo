/** @file stokes_INS_example.cpp
 * 
 * Steady Stokes example using gsIncompressibleFlow

    Author(s): M. Riegler
*/

#include <gismo.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsFlowBndEvaluators.h>

using namespace gismo;

template<class T, int MatOrder> void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt);

// Global Typedefs
typedef gsExprAssembler<>::geometryMap geometryMap;
typedef gsExprAssembler<>::variable variable;
typedef gsExprAssembler<>::space space;
typedef gsExprAssembler<>::solution solution;

int main(int argc, char *argv[])
{
  typedef gsGMRes<real_t> LinSolver;

  bool steady = true;

  std::string inputFile = "gismo/convergence_studies/moeller.xml";

  std::string path = gsFileManager::find(inputFile);
  if ( path.empty() )
  {
      gsWarn<<"Input file not found, quitting.\n";
      return 1;
  }

  int numRefine = 0;
  int numElevate = 0;

  real_t viscosity = 1.0;

  std::string matFormation = "EbE";
  std::string precond = "MSIMPLER_FdiagEqual";

  bool use_stabilization, use_pspg{false}, use_gls{false}, equal_order_bases{false};
  bool plot{true};

  gsCmdLine cmd("Solves the Navier-Stokes problem in a given domain (step, cavity, blade profile).");

  cmd.addString("f", "file", "Input XML file", inputFile);
  cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
  cmd.addInt("e", "degElevate", "Number of degree elevations (performed before h-refinement)", numElevate);

  cmd.addSwitch("pspg", "Using PSPG stabilization", use_pspg);
  cmd.addSwitch("gls", "Use GLS stabilization", use_gls);
  cmd.addSwitch("equal-order", "Use equal order bases for velocity and pressure", equal_order_bases);
  cmd.addSwitch("no-plot", "Do not plot the solution", plot);

  try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

  use_stabilization = use_pspg or use_gls;

  gsFileData<> fd(inputFile);
  gsMultiPatch<> patches;
  gsBoundaryConditions<> bcInfo, pBcInfo;
  gsFunctionExpr<> f; // external force

  fd.getId(0, patches);   // id=0: multipatch domain
  fd.getId(1, bcInfo);    // id=1: all boundary conditions
  // fd.getId(2, pBcInfo);     // id=2: pressure boundary conditions
  fd.getId(100, f);         // id=100: source function

  gsInfo << "Solving stationary Stokes.\n";
  gsInfo << patches;
  gsInfo << "viscosity = " << viscosity << "\n";
  gsInfo << "source function = " << f << "\n";

  // ========================================= Define basis ========================================= 

  // Define discretization space by refining the basis of the geometry
  gsMultiBasis<> basis(patches);
  basis.degreeElevate(numElevate);

  for (int r = 0; r < numRefine; ++r) {
    basis.uniformRefine();
  }

  std::vector< gsMultiBasis<> >  discreteBases;
  discreteBases.push_back(basis); // basis for velocity
  discreteBases.push_back(basis); // basis for pressure
  if (!equal_order_bases) {
    discreteBases[0].degreeElevate(1); // elevate the velocity space (Taylor-Hood element type)
  }
  
  gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
  gsFlowSolverParams<real_t> params(NSpde, discreteBases);
  params.options().setSwitch("quiet", false);
  params.options().setString("assemb.loop", matFormation);

  gsOptionList solveOpt;
  solveOpt.addInt("geo", "", 0);
  // solveOpt.addInt("maxIt", "", 600);
  solveOpt.addInt("plotPts", "", 10000);
  // solveOpt.addInt("animStep", "", animStep);
  // solveOpt.addReal("tol", "", tol);
  solveOpt.addSwitch("plot", "", plot);
  // solveOpt.addSwitch("animation", "", animation);
  solveOpt.addSwitch("plotMesh", "", true);
  // solveOpt.addSwitch("stokesInit", "", stokesInit);
  solveOpt.addString("id", "", "");

  // Steady without any iterations
  solveOpt.setString("id", "steady");
  params.options().setString("lin.solver", "direct");

  gsINSSolverSteady<real_t, ColMajor> NSsolver(params);

  gsInfo << "\n----------\n";
  gsInfo << "Solving the steady problem with direct linear solver.\n";

  solveProblem(NSsolver, solveOpt);

  // ------------------------------------------------------------------------
  // ___ ___ ___  __  ___    ___ __  _   ____  _ _    __ _____ _  __  __  _  
  // | __| _ \ _ \/__\| _ \  / _//  \| | / _/ || | |  /  \_   _| |/__\|  \| | 
  // | _|| v / v / \/ | v / | \_| /\ | || \_| \/ | |_| /\ || | | | \/ | | ' | 
  // |___|_|_\_|_\\__/|_|_\  \__/_||_|___\__/\__/|___|_||_||_| |_|\__/|_|\__| 
  // ------------------------------------------------------------------------
  gsExprAssembler<> expr_assembler(2, 2);
  expr_assembler.setIntegrationElements(discreteBases[0]);
  geometryMap geoMap = expr_assembler.getMap(patches);
  bcInfo.setGeoMap(patches);
  gsExprEvaluator<> expression_evaluator(expr_assembler);

  // Error calculation
  gsField<> velocity_field = NSsolver.constructSolution(0);
  gsField<> pressure_field = NSsolver.constructSolution(1);
  auto solution_vector = NSsolver.getSolution();


  // Trial to compute the error
  gsFunctionExpr<> vel_analytical_func, p_analytical_func;
  fd.getId(101, vel_analytical_func);
  fd.getId(102, p_analytical_func);
  // Type: gsExprHelper<T>::variable
  variable velocity_analytical = expression_evaluator.getVariable(vel_analytical_func);
  variable pressure_analytical = expression_evaluator.getVariable(p_analytical_func);
  auto velocity_error = velocity_field.distanceL2(velocity_analytical.source());
  auto pressure_error = pressure_field.distanceL2(pressure_analytical.source());
  // auto velocity_error = velocity_field.distanceL2(velocity_analytical, discreteBases[0]);   // possible additional parameters: isFunc_param=false or gsMultitBasis B plus isFunc_param



  // real_t velocity_error, pressure_error;
  // velocity_error = math::sqrt(expression_evaluator.integral(
  //   (vel_analytical - velocity_field).sqNorm() * meas(geoMap)
  // ));
  // pressure_error = math::sqrt(expression_evaluator.integral(
  //   (p_analytical - pressure_field).sqNorm() * meas(geoMap)
  // ));

  gsInfo << "Errors:\nVelocity: " << velocity_error << '\n'
        << "Pressure: " << pressure_error << '\n';


  return EXIT_SUCCESS;
}

template<class T, int MatOrder>
void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt)
{
    gsStopwatch clock;

    // ------------------------------------
    // prepare strings for output filenames

    bool plot = opt.getSwitch("plot");
    std::string geoStr = "customGeo";
    std::string id = opt.getString("id");
    if (plot)
    {
        index_t dim = NSsolver.getParams()->getPde().domain().geoDim();
        std::string dimStr = util::to_string(dim) + "D";
    }

    // ------------------------------------
    // solve problem

    gsInfo << "\ninitialization...\n";
    NSsolver.initialize();

    gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";
    
    NSsolver.solveStokes(); 
    // NSsolver.solve(1000, 1e-8, 0);

    real_t totalT = clock.stop();

    gsInfo << "\nAssembly time:" << NSsolver.getAssemblyTime() << "\n";
    gsInfo << "Solve time:" << NSsolver.getSolveTime() << "\n";
    gsInfo << "Solver setup time:" << NSsolver.getSolverSetupTime() << "\n";
    gsInfo << "Total solveProblem time:" << totalT << "\n\n";

    // ------------------------------------
    // plot

    if (plot) 
    {
        gsField<> velocity = NSsolver.constructSolution(0);
        gsField<> pressure = NSsolver.constructSolution(1);

        int plotPts = opt.getInt("plotPts");

        gsInfo << "Plotting in Paraview...";
        gsWriteParaview<>(velocity, geoStr + "_" + "_velocity", plotPts, opt.getSwitch("plotMesh"));
        gsWriteParaview<>(pressure, geoStr + "_" + "_pressure", plotPts);
        // plotQuantityFromSolution("divergence", velocity, geoStr + "_" + id + "_velocityDivergence", plotPts);
        gsInfo << " done.\n";
    }
}