/** @file fitting_example.cpp

    @brief Demonstrates adaptive fitting of data samples

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Options with default values
    index_t numURef   = 0;  // r
    index_t numKnots   = 0; // n
    index_t nx   = 5; // a
    index_t ny   = 2; // b
    index_t deg_x     = 2;  // x
    index_t deg_y     = 2;  // y
    real_t lambda = 1e-06;  // s
    std::string fn = "../filedata/fitting/shipHullPts55_scale01.xml"; // d, string file input data

    // Reading options from the command line
    gsCmdLine cmd("Fit parametrized sample data with a surface patch. Expected input file is an XML "
            "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
            "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
            "3 x N matrix. Every column represents a point (x,y,z) in space.");

    cmd.addInt("x", "deg_x", "degree in x direction", deg_x);
    cmd.addInt("y", "deg_y", "degree in y direction", deg_y);
    cmd.addReal("s", "lambda", "smoothing coefficient", lambda);
    
    cmd.addInt("r", "urefine", "initial uniform refinement steps", numURef);
    cmd.addInt("n", "iknots", "number of interior knots in each direction", numKnots);
    cmd.addInt("a", "uknots", "number of interior knots in u-direction", nx);
    cmd.addInt("b", "vknots", "number of interior knots in v-direction", ny);
    
    cmd.addString("d", "data", "Input sample data", fn);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (deg_x < 1)
    { gsInfo << "Degree x must be positive.\n";  return 0;}
    if (deg_y < 1)
    { gsInfo << "Degree y must be positive.\n"; return 0;}
    

    //! [Read data]
    // Surface fitting
    // Expected input is a file with matrices with:
    // id 0:  u,v   -- parametric coordinates, size 2 x N
    // id 1:  x,y,z -- corresponding mapped values, size 3 x N
    gsFileData<> fd_in(fn);
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);
    //! [Read data]

    gsWriteParaviewPoints(uv, "uv");
    gsWriteParaviewPoints(xyz, "xyz");

    // This is for outputing an XML file, if requested
    gsFileData<> fd;

    // Check if matrix sizes are OK
    GISMO_ASSERT( uv.cols() == xyz.cols() && uv.rows() == 2 && xyz.rows() == 3, "Wrong input");

    // Determine the parameter domain by mi/max of parameter values
    real_t u_min = uv.row(0).minCoeff(),
        u_max = uv.row(0).maxCoeff(),
        v_min = uv.row(1).minCoeff(),
        v_max = uv.row(1).maxCoeff();

    // Create knot-vectors without interior knots
    if( nx < 0)
      nx = numKnots;
    if( ny < 0)
      ny = numKnots;
    gsKnotVector<> u_knots (u_min, u_max, nx, deg_x+1 ) ;
    gsKnotVector<> v_knots (v_min, v_max, ny, deg_y+1 ) ;

    // Create a tensor-basis nad apply initial uniform refinement
    gsTensorBSplineBasis<2> tbasis( u_knots, v_knots );
    tbasis.uniformRefine( (1<<numURef)-1 );

    // Print settings summary
    gsInfo<<"Fitting "<< xyz.cols() <<" samples.\n";
    gsInfo<<"----------------\n";
    gsInfo<<"Bi-degree: ("<< deg_x<< ", " << deg_y << ")\n";
    gsInfo<<"Knots in x: "<< nx<<".\n";
    gsInfo<<"Knots in y: "<< ny<<".\n";
    gsInfo<<"Smoothing parameter: "<< lambda<<".\n";
    gsInfo<<"----------------\n";


    // Create fitting object
    gsFitting<real_t> fitter( uv, xyz, tbasis);
    fitter.compute(lambda);
    fitter.computeErrors();

    gsInfo << "------------------------------------\n";
    gsInfo << "computeErrors():\n";
    gsInfo << "Max error: " << fitter.maxPointError() << "\n";
    gsInfo << "Min error: " << fitter.minPointError() << "\n";
    const std::vector<real_t> & errors = fitter.pointWiseErrors();
    gsInfo << "First point-wise error: " << errors[0] << "\n";
    
    
    gsInfo << "------------------------------------\n";
    gsInfo << "pointWiseErrors(parameters,points):\n";
    gsMatrix<> ptwise = fitter.pointWiseErrors(fitter.returnParamValues(), xyz);
    gsInfo << "First point-wise error: " << ptwise(0,0) << "\n";
    gsInfo << "comparison:" << ptwise(0,0) - errors[0] << "\n";


    gsInfo << "------------------------------------\n";
    gsInfo << "computeErrors(parameters,points):\n";
    std::vector<real_t> min_max_mse = fitter.computeErrors(fitter.returnParamValues(), xyz);
    gsInfo << "Max error: " << min_max_mse[1]<< "\n";
    gsInfo << "Min error: " << min_max_mse[0] << "\n";
    gsInfo << "MSE error: " << min_max_mse[2] << "\n";

    gsInfo << "comparison of MAX error: " << fitter.maxPointError() - min_max_mse[1] << "\n";
    gsInfo << "comparison of MIN error: " << fitter.minPointError() - min_max_mse[0] << "\n";
    
    gsInfo << "------------------------------------\n";
    // Create a second fitting object for comparison
    gsFitting<real_t> fitter2( uv, xyz, tbasis);
    fitter2.compute(lambda);
    fitter2.computeMaxNormErrors();
    gsInfo << "computeMaxNormErrors():\n";
    gsInfo << "Max error: " << fitter2.maxPointError() << "\n";
    gsInfo << "Min error: " << fitter2.minPointError() << "\n";
    const std::vector<real_t> & errors2 = fitter2.pointWiseErrors();
    gsInfo << "First point-wise error: " << errors2[0] << "\n";

    return 0;
}
