/** @file CDR_example.cpp
 * 
 * Short convection example to test the convection of the stokes_temperature_example
 * against. The results should be the same
*/

//! [Include namespace]
# include <gismo.h>

using namespace std;
using namespace gismo;
//! [Include namespace]

template<class T>
void printMatrix(gsMatrix<T> mat) {
   gsDebug << "---------------------------------\n";
  for (size_t i = 0; i < mat.rows(); i++) {
    gsDebug << mat.at(i) << '\n';
  }
  gsDebug << "---------------------------------\n";
}

int main(int argc, char *argv[])
{
   //! [Parse command line]
   bool plot = true;
   std::string inputFile = "pde/poiseuille_temperature.xml";
   int numRefine = 0;
   int numElevate = 0;


   gsCmdLine cmd("Example for solving a convection-diffusion problem.");
   cmd.addString("f", "file", "Input XML file", inputFile);
   cmd.addSwitch("no-plot", "Create a ParaView visualization file with the solution", plot);
   cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
   cmd.addInt("e", "degElevate", "Number of degree elevations (performed before h-refinement)", numElevate);
   try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
   //! [Parse command line]

   // --------------- specify exact solution and right-hand-side ---------------

   //! [Function data]
   // Define source function
   gsFunctionExpr<> rhs("0",2);

   // diffusion coefficient:
   gsFunctionExpr<> coeff_diff("0.0","0.0","0.0","0.0",2);
   // convection coefficient:
   gsFunctionExpr<> coeff_conv("y*(1-y)","0.0",2);
   // reaction coefficient:
   gsFunctionExpr<> coeff_reac("0",2);
   //! [Function data]

   // Print out source function and solution
   gsInfo<<"Source function " << rhs << "\n";


   // --------------- read geometry from file ---------------

   gsFileData<> fd(inputFile);
   gsMultiPatch<> patches;
   gsBoundaryConditions<> bcInfo;
   gsFunctionExpr<> f; // external force

   fd.getId(0, patches);   // id=0: multipatch domain
   fd.getId(66, bcInfo);    // id=66: temperature boundary conditions
   fd.getId(100, f);         // id=100: source function

   //! [GetGeometryData]
   gsInfo << "The domain is a "<< patches <<"\n";

   // --------------- define Pde ---------------
   //! [definePde]
   gsConvDiffRePde<real_t> cdrPde(patches, bcInfo, & coeff_diff,& coeff_conv, & coeff_reac, & rhs);
   //! [definePde]


   // --------------- set up basis ---------------

   //! [GetBasisFromTHB]
   // Copy basis from the geometry
   gsMultiBasis<> bases( patches );
   //! [GetBasisFromTHB]


   //! [Refinements]
   bases.degreeElevate(numElevate);
   for (int r = 0; r < numRefine; ++r) {
      bases.uniformRefine();
   }

   //! [constructAssembler]
   // Construct assembler
   gsCDRAssembler<real_t> cdrAss( cdrPde, bases);
   // Set stabilization flag to 1 = SUPG
   cdrAss.options().setInt("Stabilization", stabilizerCDR::SUPG);
   // Compute Dirichlet values by L2-projection
   // Caution: Interpolation does not work for locally refined (T)HB-splines!
   cdrAss.options().setInt("DirichletValues",dirichlet::l2Projection);
   //! [constructAssembler]


   //! [beginRefLoop]
   // --------------- solving ---------------

   //! [solverPart]
   // Generate system matrix and load vector
   cdrAss.assemble();

   auto matValues = cdrAss.matrix().valuePtr();
   auto outerIndex = cdrAss.matrix().outerIndexPtr();
   auto innerIndex = cdrAss.matrix().innerIndexPtr();

   // Solve the system
   gsMatrix<real_t> solVector =
      gsSparseSolver<>::BiCGSTABILUT( cdrAss.matrix() ).solve( cdrAss.rhs() );

   // Construct the solution as a scalar field
   gsField<> solField;
   solField = cdrAss.constructSolution(solVector);
   //! [solverPart]

   //! [Plot in Paraview]
   if( plot )
   {
       // Run paraview
       gsWriteParaview<>(solField, "Temperature_convection.pvd", 1000, false);
      //  gsFileManager::open("Temperature_convection.pvd");
   }
   //! [Plot in Paraview]
   else
   {
       gsInfo<<"Done. No output created, re-run with --plot to get a ParaView "
               "file containing Plotting image data.\n";
   }
   return EXIT_SUCCESS;

}// end main
