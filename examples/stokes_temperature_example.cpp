/** @file stokes_temperature_example.cpp
 * 
 * Steady Stokes example using gsIncompressibleFlow and afterwards advecting temperature
 * field with the fluid flow (including a source term stemming from viscous dissipation
 * in shear-thinning fluids)

    Author(s): M. Riegler
*/

#include <gismo.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsFlowBndEvaluators.h>

using namespace gismo;

template<class T> void solveHeatProblem(gsMultiPatch<T> patches, gsField<T> velocityField, index_t dim,
                                        gsBoundaryConditions<T> temperatureBcInfo, int numRefine,
                                        int numElevate);
template<class T, int MatOrder> void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt);

template<class T>
void printMatrix(gsMatrix<T> mat) {
   gsDebug << "---------------------------------\n";
  for (size_t i = 0; i < mat.rows(); i++) {
    gsDebug << mat.at(i) << '\n';
  }
  gsDebug << "---------------------------------\n";
}

// Global Typedefs
typedef gsExprAssembler<>::geometryMap geometryMap;
typedef gsExprAssembler<>::variable variable;
typedef gsExprAssembler<>::space space;
typedef gsExprAssembler<>::solution solution;

int main(int argc, char *argv[])
{
  typedef gsGMRes<real_t> LinSolver;

  bool steady = true;

  std::string inputFile = "microstructures/microstructure.xml";

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
  solveOpt.addSwitch("plotMesh", "", false);
  // solveOpt.addSwitch("stokesInit", "", stokesInit);
  solveOpt.addString("id", "", "");

  // Steady without any iterations
  solveOpt.setString("id", "steady");
  params.options().setString("lin.solver", "direct");

  gsINSSolverSteady<real_t, ColMajor> NSsolver(params);

  gsInfo << "\n----------\n";
  gsInfo << "Solving the steady problem with direct linear solver.\n";

  solveProblem(NSsolver, solveOpt);

// --------------------------------------------------------------------
//    _  _  ____   __  ____    ____  ____   __  ____  __    ____  _  _ 
//   / )( \(  __) / _\(_  _)  (  _ \(  _ \ /  \(  _ \(  )  (  __)( \/ )
//   ) __ ( ) _) /    \ )(     ) __/ )   /(  O )) _ (/ (_/\ ) _) / \/ \
//   \_)(_/(____)\_/\_/(__)   (__)  (__\_) \__/(____/\____/(____)\_)(_/
// --------------------------------------------------------------------
  gsField<> velocityField = NSsolver.constructSolution(0);

  index_t dim = NSsolver.getParams()->getPde().domain().geoDim();

  // Get boundary conditions
  gsBoundaryConditions<> temperatureBcInfo;
  fd.getId(66, temperatureBcInfo);

  solveHeatProblem<>(patches, velocityField, dim, temperatureBcInfo, numRefine, numElevate);

  return EXIT_SUCCESS;
}

template<class T>
void solveHeatProblem(gsMultiPatch<T> patches, gsField<T> velocityField, index_t dim,
                      gsBoundaryConditions<T> temperatureBcInfo, int numRefine,
                      int numElevate) {
  const gsFunctionSet<>& velocityFunction = velocityField.fields();
  // gsFunctionExpr<> coeff_diffusion(diff_coeff_str, dim);
  gsFunctionExpr<> coeff_diffusion("0.0","0","0","0.0",2);
  // Reaction term
  gsFunctionExpr<> coeff_reaction("0.0", dim);
  // TODO; for now no rhs term, but velocity-dependent viscous dissipation as source term
  gsFunctionExpr<> rhs("0.0", dim);
  
  
  

  // Define PDE and assembler
  gsConvDiffRePde<> cdrPde(patches, temperatureBcInfo, &coeff_diffusion,
      &velocityFunction, &coeff_reaction, &rhs);
  gsMultiBasis<> functionBasisTemperature(patches);
  functionBasisTemperature.setDegree(functionBasisTemperature.maxCwiseDegree() + numElevate);
  for (int r = 0; r < numRefine; r++) {
      functionBasisTemperature.uniformRefine();
  }
  gsCDRAssembler<> cdrAss(cdrPde, functionBasisTemperature);
  cdrAss.options().setInt("Stabilization", stabilizerCDR::SUPG);
  cdrAss.options().setInt("DirichletValues", dirichlet::l2Projection);

  cdrAss.assemble();

  gsMatrix<> temperatureSolVector = gsSparseSolver<>::BiCGSTABILUT(cdrAss.matrix()).solve(cdrAss.rhs());

  gsField<> temperatureField = cdrAss.constructSolution(temperatureSolVector);

  gsInfo << "After having solved the heat problem\n";

  // Plotting
  gsWriteParaview<>(temperatureField, "customGeo_temperature", 1000, false);

  gsInfo << "Finished exporting the temperature\n";
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

        // gsParaviewCollection collection("ParaviewOutput/" + geoStr, )

        gsInfo << " done.\n";
    }
}