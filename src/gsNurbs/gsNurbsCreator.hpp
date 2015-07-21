/** @file gsNurbsCreator.hpp

    @brief Provides implementation of the NurbsCreator struct.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsMultiPatch.h>

#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorNurbs.h>
#include <gsNurbs/gsNurbs.h>
#include <gsNurbs/gsBSpline.h>

namespace gismo
{

/*
   @brief Class gsNurbsCreator provides some simple examples of Nurbs Geometries
*/

template<class T>
gsTensorBSpline<3,T> * gsNurbsCreator<T>::lift3D( gsTensorBSpline<2,T> const & geo, T z)
{
    gsKnotVector<T> KV(0, 1, 0, 2);
    const int sz = geo.basis().size();

    gsMatrix<T> newcoefs( 2*sz, geo.geoDim() ) ;

    // Copy coefficients
    newcoefs.topRows(sz)    =
        newcoefs.bottomRows(sz) = geo.coefs();

    // Embed in 3D if needed
    if (newcoefs.cols() == 2 )
    {
        newcoefs.conservativeResize( Eigen::NoChange, 3);
        newcoefs.col(2).setZero();
    }

    // Lift
    newcoefs.col(2).bottomRows(sz).array() += z;

    return new gsTensorBSpline<3,T>(geo.basis().knots(0),
                                    geo.basis().knots(1),
                                    KV, give(newcoefs) );
}

    
template<class T> gsTensorBSpline<4,T> * gsNurbsCreator<T>::lift4D( gsTensorBSpline<3,T> const & geo, T z)
{
    gsKnotVector<T> KV(0, 1, 0, 2);
    const int sz = geo.basis().size();

    gsMatrix<T> newcoefs( 3*sz, geo.geoDim() ) ;

    // Copy coefficients
    newcoefs.topRows(sz)    =
        newcoefs.bottomRows(sz) = geo.coefs();

    // Embed in 4D if needed
    if (newcoefs.cols() == 3 )
    {
        newcoefs.conservativeResize( Eigen::NoChange, 4);
        newcoefs.col(3).setZero();
    }

    // Lift
    newcoefs.col(3).bottomRows(sz).array() += z;

    return new gsTensorBSpline<4,T>(geo.basis().knots(0),
                                    geo.basis().knots(1),
                                    geo.basis().knots(2),
                                    KV, give(newcoefs) );
}

/* * Computes a set of control points, weights, and knots that define an order-3 circular arc centered at the origin
    \param X Defines the X axis of the plane containing the arc
    \param Y Defines the Y axis of the plane containing the arc
    \param StartAngle Start angle of the arc in radians
    \param EndAngle End angle of the arc in radians
    \param Segments The number of NURBS segments in the resulting arc
    \param Knots Output container for the resulting arc knot vector
    \param Weights Output container for the resulting arc control point weights
    \param ControlPoints Output container for the resulting arc control point positions

    template<class T> gsNurbsCreator<T> * gsNurbsCreator<T>::circularArc(const vector3& X, const vector3& Y,
    const T StartAngle, const T EndAngle,
    const unsinged Segments = 1)
    {
    gsKnotVector<T> Knots;
    gsMatrix<T>     Weights;
    gsMatrix<T>     ControlPoints;

    const T theta = (EndAngle - StartAngle) / static_cast<T>(Segments);
    const T weight = std::cos(std::fabs(theta) * 0.5);

    Knots.clear();
    Knots.insert(Knots.end(), 3, 0);
    for(uint_t i = 1; i != Segments; ++i)
    Knots.insert(Knots.end(), 2, i);
    Knots.insert(Knots.end(), 3, Segments);

    Weights.clear();
    Weights.push_back(1.0);
    for(uint_t i = 0; i != Segments; ++i)
    {
    Weights.push_back(weight);
    Weights.push_back(1);
    }

    ControlPoints.clear();
    ControlPoints.push_back(k3d::to_point(std::cos(StartAngle) * X + std::sin(StartAngle) * Y));
    for(uint_t i = 0; i != Segments; ++i)
    {
    const T a0 = math::mix(StartAngle, EndAngle, static_cast<T>(i) / static_cast<T>(Segments));
    const T a2 = math::mix(StartAngle, EndAngle, static_cast<T>(i+1) / static_cast<T>(Segments));

    const point3 p0(k3d::to_point(std::cos(a0) * X + std::sin(a0) * Y));
    const point3 p2(k3d::to_point(std::cos(a2) * X + std::sin(a2) * Y));

    const point3 t0(k3d::to_point(-std::sin(a0) * X + std::cos(a0) * Y));
    const point3 t2(k3d::to_point(-std::sin(a2) * X + std::cos(a2) * Y));

    point3 p1;
    intersect_lines(p0, to_vector(t0), p2, to_vector(t2), p1);

    ControlPoints.push_back(p1);
    ControlPoints.push_back(p2);
    }
    }
*/

template<class T> gsBSpline<T> * gsNurbsCreator<T>::BSplineUnitInterval(int deg)
{
    gsKnotVector<T> KV(0,1, 0,deg+1);
    gsMatrix<T> C(deg+1, 1);
    for (int i = 0; i <= deg; ++i)
        C(i) = i / static_cast<T>(deg);
    return new gsBSpline<T>(KV, give(C));
}

/// 2d-rectange [low_x,upp_x] x [low_y,upp_y], rotated by \a turndeg degrees.
template<class T> gsTensorBSpline<2,T> * 
gsNurbsCreator<T>::BSplineRectangle( T const & low_x,
                                     T const & low_y,
                                     T const & upp_x,
                                     T const & upp_y, T const & turndeg)
{

    gsKnotVector<T> KV (0,1,0,3);
    gsMatrix<T> C(4,2);

    const T pi = 3.1415926535897932384626433832795;

    T r = turndeg / 180 * pi;

    C <<  low_x, low_y ,
        upp_x, low_y,
        low_x, upp_y,
        upp_x, upp_y;

    T tx;
    T ty;
    for(int i =0; i < 4; i++)
    {
        tx = C(i,0); ty = C(i,1);
        C(i,0) = std::cos(r) * tx - std::sin(r) * ty;
        C(i,1) = std::sin(r) * tx + std::cos(r) * ty;
    }

    gsMatrix<T> D(9,2);
    D.setZero();
    D(0,0) = C(0,0); D(0,1) = C(0,1);
    D(2,0) = C(1,0); D(2,1) = C(1,1);
    D(6,0) = C(2,0); D(6,1) = C(2,1);
    D(8,0) = C(3,0); D(8,1) = C(3,1);

    D(1,0) = (C(0,0)+C(1,0))/2; D(1,1) = (C(0,1)+C(1,1))/2;
    D(3,0) = (C(0,0)+C(2,0))/2; D(3,1) = (C(0,1)+C(2,1))/2;
    D(5,0) = (C(3,0)+C(1,0))/2; D(5,1) = (C(3,1)+C(1,1))/2;
    D(7,0) = (C(2,0)+C(3,0))/2; D(7,1) = (C(2,1)+C(3,1))/2;
    D(4,0) = (C(0,0)+C(3,0))/2; D(4,1) = (C(0,1)+C(3,1))/2;

    return new gsTensorBSpline<2,T>(KV,KV, give(D));

}


template<class T> gsTensorBSpline<2,T> * 
gsNurbsCreator<T>::BSplineRectangleWithPara( T low_x, T low_y, T upp_x, T upp_y)
{
    gsKnotVector<T> KVx (low_x, upp_x, 0, 2);
    gsKnotVector<T> KVy (low_y, upp_y, 0, 2);

    gsMatrix<T> C(4,2);

    C << low_x, low_y,
         upp_x, low_y,
         low_x, upp_y,
         upp_x, upp_y;

    return new gsTensorBSpline<2,T>(KVx, KVy, give(C));
}


/// Square of side \a r, with lower left corner at (x,y)
template<class T> gsTensorBSpline<2,T> * 
gsNurbsCreator<T>::BSplineSquare( T const & r, 
                                  T const & x,
                                  T const & y)
{
    gsKnotVector<T> KV (0,1,0,2) ;
    gsMatrix<T> C(4,2) ;

    C <<  0 , 0 , 1 , 0
        , 0 , 1 , 1 , 1 ;
    // scale
    C *= r;

    // translate
    C.col(0).array() += x;
    C.col(1).array() += y;

    return new gsTensorBSpline<2,T>(KV,KV, give(C));
};
/// Creates a \em n X \em m rectangle multipatch consisting of B-splines squares
/// with lower left corner at at (lx,ly).
/// The numbering of the patches are (for \em n = 4 and \em m = 3):
/// --|---|---|--
/// --|---|---|--
/// 2 | 5 | 8 |11
/// 1 | 4 | 7 |10
/// 0 | 3 | 6 |9
/// \param n number of squares in x-direction.
/// \param m number of squares in y-direction.
/// \param r with length of the side of the squares.
/// \param lx x-coordinate for lower left corner of the rectangle.
/// \param ly y-coordinate for lower left corner of the rectangle.
template<class T> gsMultiPatch<T> * 
gsNurbsCreator<T>::BSplineSquareGrid(int n, int m, 
                                     T const & r,
                                     T const & lx, 
                                     T const & ly)
{
    gsMultiPatch<T> * mp = new gsMultiPatch<T>;

    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
        {
            mp->addPatch(BSplineSquare(r,lx + r*i ,ly + r*j)) ;
        }
    mp->computeTopology();
    return mp;
};


template<class T> gsTensorBSpline<2,T> * 
gsNurbsCreator<T>::BSplineSquare( gsMatrix<T> const & Box)
{
    gsKnotVector<T> KV (0,1,0,2) ;
    gsMatrix<T> C(4,2) ;
    C << Box(0,0) , Box(1,0), Box(0,1),Box(1,0),
        Box(0,0) , Box(1,1), Box(0,1),Box(1,1)  ;

    return new gsTensorBSpline<2,T>(KV,KV, give(C));
};


// Note: this can probably be removed once we have degree elevation for tensor B-splines.
//
/// The unit square represented as a tensor B-spline of degree \a deg
template<class T> gsTensorBSpline<2,T> * 
gsNurbsCreator<T>::BSplineSquare(int deg)
{
    const int n = (deg + 1) * (deg + 1);        // number of basis functions

    gsMatrix<T> C(n, 2);

    index_t r = 0;

    for (int j = 0; j <= deg; ++j)
        for (int i = 0; i <= deg; ++i)
        {
            C(r, 0) = ((T) i) / deg;
            C(r, 1) = ((T) j) / deg;
            ++r;
        }

    gsKnotVector<T> KV(0,1, 0, deg+1);
    return new gsTensorBSpline<2,T>(KV,KV, give(C));
}


template<class T> gsTensorBSpline<3,T> * 
gsNurbsCreator<T>::BSplineCube( T const & r, T const & x,
                                T const & y, T const & z)
{
    gsKnotVector<T> KV (0,1,0,2) ;
    gsMatrix<T> C(8,3) ;
    C <<  -0.5 , -0.5 , -0.5, 0.5 , -0.5 , -0.5
        , -0.5 , 0.5 , -0.5, 0.5 , 0.5 , -0.5
        , -0.5 , -0.5 , 0.5, 0.5 , -0.5 , 0.5
        , -0.5 , 0.5 , 0.5, 0.5 , 0.5 , 0.5 ;
    C *= r;
    C.col(0).array() += x;
    C.col(1).array() += y;
    C.col(2).array() += z;

    return new gsTensorBSpline<3,T>(KV,KV,KV, give(C));
};


// Note: this can probably be removed once we have degree elevation for tensor B-splines.
//
/// The unit cube represented as a tensor B-spline of degree \a deg
template<class T> gsTensorBSpline<3,T> * 
gsNurbsCreator<T>::BSplineCube(int deg)
{
    const int n = (deg + 1) * (deg + 1) * (deg + 1);        // number of basis functions

    gsMatrix<T> C(n, 3);

    index_t r = 0;

    for (int k = 0; k <= deg; ++k)
        for (int j = 0; j <= deg; ++j)
            for (int i = 0; i <= deg; ++i)
            {
                C(r, 0) = ((T) i) / deg;
                C(r, 1) = ((T) j) / deg;
                C(r, 2) = ((T) k) / deg;
                ++r;
            }

    gsKnotVector<T> KV(0,1, 0, deg+1);
    return new gsTensorBSpline<3,T>(KV,KV,KV, give(C));
}


template<class T> gsTensorBSpline<3,T> * 
gsNurbsCreator<T>::BSplineHalfCube( T const & r, T const & x,
                                    T const & y, T const & z)
{
    gsKnotVector<T> KV (0,1,0,2) ;
    gsMatrix<T> C(8,3) ;
    C <<  -0.5 , -0.5 , -0.5, 0.5 , -0.5 , -0.5
        , -0.5 , 0.5  , -0.5, 0   , 0    , -0.5
        , -0.5 , -0.5 , 0.5, 0.5  , -0.5 , 0.5
        , -0.5 , 0.5  , 0.5, 0    , 0    , 0.5 ;
    C *= r;
    C.col(0).array() += x;
    C.col(1).array() += y;
    C.col(2).array() += z;

    return new gsTensorBSpline<3,T>(KV,KV,KV, give(C));
};


template<class T> gsTensorNurbs<3,T> * 
gsNurbsCreator<T>::NurbsCube( T const & r, T const & x,
                              T const & y, T const & z)
{
    gsKnotVector<T> KV (0,1,0,2) ;
    gsMatrix<T> C(8,3) ;
    C <<  0 , 0 , 0, 1 , 0 , 0  // Cube
        , 0 , 1 , 0, 1 , 1 , 0
        , 0 , 0 , 1, 1 , 0 , 1
        , 0 , 1 , 1, 1 , 1 , 1 ;
        
    C.col(0).array() += x;
    C.col(1).array() += y;
    C.col(2).array() += z;        

    return new gsTensorNurbs<3,T>(KV,KV,KV, give(C));
};

template<class T> gsTensorNurbs<2,T> * 
gsNurbsCreator<T>::NurbsQuarterAnnulus( T const & r0, T const & r1)
{
    gsKnotVector<T> KVx (0,1,0,2) ;
    gsKnotVector<T> KVy (0,1,0,3) ;
    gsMatrix<T> C(6,2) ;
    C <<  r0 , 0  ,  r1, 0
        , r0 , r0 ,  r1, r1
        , 0 ,  r0 ,  0 , r1 ;

    // Set weights
    gsMatrix<T> ww(6,1) ;
    ww.setOnes();
    ww.at(2)= 0.707106781186548 ;
    ww.at(3)= 0.707106781186548 ;

    return new gsTensorNurbs<2,T>(KVx,KVy, give(C), give(ww));
}

template<class T> gsTensorNurbs<3,T> * 
gsNurbsCreator<T>::BSplineSaddle()
{
    return NULL;
}

/// Inexact annulus using B-splines
template<class T> gsGeometry<T> * 
gsNurbsCreator<T>::BSplineQuarterAnnulus(int const & deg)
{
    gsGeometry<T> * quann = gsNurbsCreator<T>::NurbsQuarterAnnulus();

    gsKnotVector<T> KV1(0,1,     0,     2);
    gsKnotVector<T> KV2(0,1, deg-2, deg+1);

    gsTensorBSplineBasis<2,T> tbsp (new gsBSplineBasis<T>(KV1), new gsBSplineBasis<T>(KV2));
    const gsMatrix<T> pts = tbsp.anchors();
    gsGeometry<T> * approxGeom = tbsp.interpolateData(quann->eval(pts),pts);
    delete quann;

    return approxGeom;
}

/// Fat annulus using B-splines, discarding the weights of the exact NURBS
/// Analytical formulation (when r0 = 1 and r1 = 2):
/// (x, y) = (1 + s - s*t*t - t*t, 2*s*t -s*t*t + 2*t - t*t)
template<class T> gsTensorBSpline<2,T> * 
gsNurbsCreator<T>::BSplineFatQuarterAnnulus( T const & r0, T const & r1)
{
    gsKnotVector<T> KVx (0,1,0,2) ;
    gsKnotVector<T> KVy (0,1,0,3) ;

    gsMatrix<T> C(6,2) ;
    C <<  r0 , 0  , r1, 0
        , r0 , r0 , r1, r1
        , 0 ,  r0 , 0 , r1 ;

    //gsMatrix<T> C(9,2) ;
    // C <<  r0 , 0  , (r0+r1)/2,     0     , r1, 0
    //     , r0 , r0 , (r0+r1)/2, (r0+r1)/2 , r1, r1
    //     , 0 ,  r0 ,      0   , (r0+r1)/2 , 0 , r1 ;

    //return new gsTensorBSpline<2,T>(KVy,KVy,2,2, C);
    return new gsTensorBSpline<2,T>(KVx,KVy, give(C));
};


template<class T> gsTensorNurbs<2,T> * 
gsNurbsCreator<T>::NurbsSphere( T const & r, T const & x, 
                                T const & y, T const & z)
{
    gsKnotVector<T> KV1 (0,2,1,3,2) ;
    gsKnotVector<T> KV2 (0,4,3,3,2) ;
    gsMatrix<T> C(45,3) ;
    C<<
        2, 0, 0, 2, 2, 0, 0, 2, 0,-2, 2, 0,-2, 0, 0,
        2, 0, 2, 2, 2, 2, 0, 0, 2, -2, 2, 2,-2,0, 2,
        0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0,0, 2,
        2, 0, 2, 2,-2, 2, 0,-2, 2,-2,-2, 2,-2,0, 2,
        2, 0, 0, 2,-2, 0, 0,-2, 0,-2,-2, 0,-2, 0, 0,
        2, 0,-2, 2,-2,-2, 0,-2,-2,-2,-2,-2,-2, 0,-2,
        0, 0,-2, 0, 0,-2, 0, 0,-2, 0, 0,-2, 0, 0,-2,
        2, 0,-2, 2, 2,-2, 0, 2,-2,-2, 2,-2,-2, 0,-2,
        2, 0, 0, 2, 2, 0, 0, 2, 0,-2, 2, 0,-2, 0, 0;

    gsMatrix<T> W(45,1) ;
    W <<          1, 0.707106781186548,                 1, 0.707106781186548,                  1,
        0.707106781186548,               0.5, 0.707106781186548,               0.5,  0.707106781186548,
        1, 0.707106781186548,                 1, 0.707106781186548,                   1,
        0.707106781186548,               0.5, 0.707106781186548,               0.5,   0.707106781186548,
        1, 0.707106781186548,                 1, 0.707106781186548,                   1,
        0.707106781186548,               0.5, 0.707106781186548,               0.5,   0.707106781186548,
        1, 0.707106781186548,                 1, 0.707106781186548,                   1,
        0.707106781186548,               0.5, 0.707106781186548,               0.5,   0.707106781186548,
        1, 0.707106781186548,                 1, 0.707106781186548,                   1;

    C.col(0).array() += x;
    C.col(1).array() += y;
    C.col(2).array() += z; 
    
    return new gsTensorNurbs<2,T>(KV1,KV2, give(C), give(W));
};


template<class T> gsNurbs<T> * 
gsNurbsCreator<T>::NurbsCircle( T const & r, T const & x, T const & y)
{
    gsKnotVector<T> KV2 (0,1,3,3,2) ;
    gsMatrix<T> C(9,2) ;
    C <<  1, 0,
        1, 1,
        0, 1,
        -1, 1,
        -1, 0,
        -1,-1,
        0,-1,
        1,-1,
        1, 0;
    C *= r;

    C.col(0).array() += x;
    C.col(1).array() += y;    
    
    gsMatrix<T> ww(9,1) ;
    ww<< 1, 0.707106781186548, 1, 0.707106781186548,1, 0.707106781186548,1, 0.707106781186548, 1 ;

    return new gsNurbs<T>(KV2, give(ww), give(C));
};

template<class T> gsBSpline<T> * 
gsNurbsCreator<T>::BSplineFatCircle( T const & r, T const & x, T const & y)
{
    gsKnotVector<T> KV2 (0,1,3,3,2) ;
    gsMatrix<T> C(9,2) ;
    C <<  1, 0,
        1, 1,
        0, 1,
        -1, 1,
        -1, 0,
        -1,-1,
        0,-1,
        1,-1,
        1, 0;
    C *= r;
    
    C.col(0).array() += x;
    C.col(1).array() += y;
    
    return new gsBSpline<T>(KV2, give(C));
};

template<class T> gsTensorBSpline<2,T> *
gsNurbsCreator<T>::BSplineFatDisk (T const & r, T const & x, T const & y)
{
    gsKnotVector<T> KV (0,1,0,3) ;
    gsMatrix<T>  C(9,2) ;
    C <<  0, -1 ,  1,-1 , 1, 0
        ,-1, -1 ,  0, 0 , 1 ,1
        ,-1 , 0 , -1, 1 , 0 ,1 ;

    C.col(0).array() += x;
    C.col(1).array() += y;
    
    return new gsTensorBSpline<2,T>(KV,KV, give(C));
}

template<class T> gsNurbs<T> *
gsNurbsCreator<T>::NurbsCurve1 (T const & r, T const & x, T const & y)
{
    gsKnotVector<T> KV2 (0,4,3,3,2) ;
    KV2.uniformRefine();
    gsMatrix<T> C(13,2) ;
    C << 0,   -1,
        0.5,   -1,
        1, -0.5,
        1,    0,
        1,  0.5,
        0.5,    1,
        0,    1,
        -0.5,    1,
        -1,  0.5,
        -1,    0,
        -1, -0.5,
        -0.5,   -1,
        0,   -1;
    C *= r;
    
    C.col(0).array() += x;
    C.col(1).array() += y;    

    gsMatrix<T> ww(13,1) ;
    ww<< 1, 0.853553, 0.853553, 1, 0.853553, 0.853553, 1, 0.853553,
        0.853553, 1, 0.853553, 0.853553, 1;

    gsNurbs<T> * nn = new gsNurbs<T>(KV2, give(C), give(ww));
    // std::cout<<" nurbs:\n " <<* nn << std::endl;
    // nn->uniformRefine();
    // std::cout<<" coefs:\n " <<* nn->coefs() << std::endl;
    // std::cout<<" weights:\n " <<* nn->weights() << std::endl;


    return nn;
};


template<class T> gsNurbs<T> *
gsNurbsCreator<T>::NurbsCurve2 (T const & r, T const & x, T const & y)
{
    gsKnotVector<T> kv(0,4,3,3,2);
    gsMatrix<T> C(9,2);
    C <<  0, -2 ,  0.5,-0.5 , 2, 0
        ,-2, -2 ,  0, 0 , 0.5 ,0.5
        ,-2 , 0 , -0.5, 0.5 , 0 ,-2 ;

    C *= r;

    C.col(0).array() += x;
    C.col(1).array() += y;
    
    gsMatrix<T> ww( 9, 1 ) ;
    ww(0)= 1;
    ww(1)= 0.20 ;
    ww(2)= 1;
    ww(3)= 0.707106781186548 ;
    ww(4)= 0.5 ;
    ww(5)= 0.20 ;
    ww(6)= 1;
    ww(7)= 0.707106781186548 ;
    ww(8)= 1 ;

    return new gsNurbs<T>(kv, give(C), give(ww));
};


template<class T> gsNurbs<T> *
gsNurbsCreator<T>::NurbsBean(T const & r, T const & x, T const & y)
{
    gsKnotVector<T> kv(0,1,12,3,1,2);
    gsMatrix<T> C(15,2);
    C <<  1,0,
        1,1,
        0,2,
        0,3,
        1,4,
        1.5,5.2,
        0.5,6,
        -1,5,
        -1.5,3,
        -2,1,
        -2,-1,
        -1,-2,
        0,-2,
        1,-1,
        1,0;
        
    C.col(0).array() += x;
    C.col(1).array() += y;         

    gsMatrix<T> ww( 15, 1 ) ;
    ww.setOnes();
    return new gsNurbs<T>(kv, give(C), give(ww));
};


template<class T> gsBSpline<T> *
gsNurbsCreator<T>::BSplineE (T const & r, T const & x, T const & y)
{
    gsKnotVector<T> kv(0,1,22,4,1,3);
    gsMatrix<T> C(26,2);
    C << -2,0,
        -2,1,
        -2,2,
        -1,3,
        2,3.3,
        2.5,2.8,
        3,2,
        2.5,1.5,
        1, 1.5,
        0.5, 1.2,
        1,1,
        2,1,
        3,0.5, 2.5,-0.5,
        0.5, -0.5,
        0.5,-1.5,
        1, -1.7,
        2, -1.5,
        3, -1.5,
        3.3, -2.2,
        3, -3,
        0, -3.3,
        -1, -3,
        -2, -2,
        -2, -1,
        -2,0;

    C.col(0).array() += x;
    C.col(1).array() += y;         

    return new gsBSpline<T>(kv, give(C));
};


template<class T> gsNurbs<T> *
gsNurbsCreator<T>::NurbsAmoebaFull(T const & r, T const & x, T const & y)
{
    gsKnotVector<T> kv(0,1,19,3,1,2);
    gsMatrix<T> C(22,2);
    C <<    -10,-2,
        -10,1,
        -9,2,
        -6.8,2.1,
        -7,4,
        -5,6.5,
        -3,6,
        -0.5,3,
        0.5,3.2,
        4,6,
        7,6,
        8,3,
        5,0,
        3.8,-1,
        5,-2,
        6,-6,
        2,-8,
        -2,-6,
        -2,-2,
        -6,-3,
        -10,-3,
        -10,-2;

    C.col(0).array() += x;
    C.col(1).array() += y;   
    
    gsMatrix<T> ww(22, 1 ) ;
    ww.setOnes();
    return new gsNurbs<T>(kv, give(ww), give(C));
}

template<class T> gsBSpline<T> *
gsNurbsCreator<T>::BSplineLineSegment(gsMatrix<T> const & p0, gsMatrix<T> const & p1 )
{
    gsKnotVector<T> kv(0,1,0,2,1,1);
    gsMatrix<T> C(2,2);

    C.row(0).noalias() =  p0.transpose();
    C.row(1).noalias() =  p1.transpose();
    return new gsBSpline<T>(kv, give(C));
}


/// L-Shaped domain represented as a tensor B-spline of degree 1
template<class T> gsTensorBSpline<2,T> * 
gsNurbsCreator<T>::BSplineLShape_p1(T r)
{
    // create knot vector [0,0, 0.5, 1,1]
    gsKnotVector<T> tK1(0,1,1,2,1);
    // create knot vector [0,0, 1,1]
    gsKnotVector<T> tK2(0,1,0,2);

    gsMatrix<T> C(6,2);

    C <<-1.0,1.0,
        -1.0,-1.0,
        1.0,-1.0,
        0.0,1.0,
        0.0,0.0,
        1.0,0.0;

    C *= r;

    return new gsTensorBSpline<2,T>(tK1,tK2, give(C));
}


/// L-Shaped domain represented as a tensor B-spline of degree 2
/// with C0-continuity across the diagonal.
template<class T> gsTensorBSpline<2,T> * 
gsNurbsCreator<T>::BSplineLShape_p2C0()
{
    // create knot vector [0,0,0, 0.5,0.5, 1,1,1]
    gsKnotVector<T> tK1(0,1,1,3,2);
    // create knot vector [0,0,0, 1,1,1]
    gsKnotVector<T> tK2(0,1,0,3);

    gsMatrix<T> C(15,2);

    C << -1.0,1.0,
        -1.0,0.0,
        -1.0,-1.0,
        0.0,-1.0,
        1.0,-1.0,
        -0.6,1.0,
        -0.55,0.0,
        -0.5,-0.5,
        0.0,-0.55,
        1.0,-0.6,
        0.0,1.0,
        0.0,0.5,
        0.0,0.0,
        0.5,0.0,
        1.0,0.0;

    return new gsTensorBSpline<2,T>(tK1,tK2, give(C));
}


/// L-Shaped domain represented as a tensor B-spline of degree 2
/// with C1-continuity and double control points at the corners.
template<class T> gsTensorBSpline<2,T> * 
gsNurbsCreator<T>::BSplineLShape_p2C1()
{
    // create knot vector [0,0,0, 0.25,0.5,0.75 1,1,1]
    gsKnotVector<T> tK1(0,1,3,3);
    // create knot vector [0,0,0, 0.5, 1,1,1]
    gsKnotVector<T> tK2(0,1,1,3);

    gsMatrix<T> C(24,2);

    C <<    -1.0, 1.0,
        -1.0, 0.2,
        -1.0, -1.0,
        -1.0, -1.0,
        0.2, -1.0,
        1.0, -1.0,
        -0.75, 1.0,
        -0.7, 0.35,
        -0.6, -0.3,
        -0.3, -0.6,
        0.35, -0.7,
        1.0, -0.75,
        -0.3, 1.0,
        -0.3, 0.5,
        -0.2, -0.05,
        -0.05, -0.2,
        0.5, -0.3,
        1.0, -0.3,
        0.0, 1.0,
        0.0, 0.5,
        0.0, 0.0,
        0.0, 0.0,
        0.5, 0.0,
        1.0, 0.0;

    return new gsTensorBSpline<2,T>(tK1,tK2, give(C));
}

template<class T> gsBSpline<T> *
gsNurbsCreator<T>::BSplineAmoeba(T const & r, T const & x, T const & y)
{
    gsKnotVector<T> kv(0,1,19,3,1,2);
    gsMatrix<T> C(22,2);

    C << -10,-2,
        -10,-1,
        -9,2,
        -6.8,2.1,
        -7,4,
        -5,6.5,
        -3,6,
        -0.5,3,
        0.5,3.2,
        4,6,
        7,6,
        8,3,
        5,0,
        3.8,-1,
        5,-2,
        6,-6,
        2,-8,
        -2,-6,
        -2,-2,
        -6,-3,
        -10,-3,
        -10,-2;

    C /= 2;

    C.col(0).array() += x;
    C.col(1).array() += y;    

    gsBSpline<T> *B =new gsBSpline<T>(kv, give(C));
    B->reverse();
    return B;
}

template<class T> gsBSpline<T> *
gsNurbsCreator<T>::BSplineAmoebaBig(T const & r, T const & x, T const & y)
{
    gsKnotVector<T> kv(0,1,19,3,1,2);
    gsMatrix<T> C(22,2);

    C << -10,-2,
        -10,-1,
        -9,2,
        -6.8,2.1,
        -7,4,
        -5,6.5,
        -3,6,
        -0.5,3,
        0.5,3.2,
        4,6,
        7,6,
        8,3,
        5,0,
        3.8,-1,
        5,-2,
        6,-6,
        2,-8,
        -2,-6,
        -2,-2,
        -6,-3,
        -10,-3,
        -10,-2;

    C.col(0).array() += x;
    C.col(1).array() += y;        

    gsBSpline<T> *B =new gsBSpline<T>(kv, give(C));
    B->reverse();
    return B;
}

template<class T> gsBSpline<T> *
gsNurbsCreator<T>::BSplineAustria(T const & r, T const & x, T const & y)
{
    gsKnotVector<T> kv(0,1,31,3,1,2);
    gsMatrix<T> C(34,2);

    C << 0.3, 2.5,
        0.7, 2.7,
        1 , 3,
        1.5, 3,
        2, 2.5,
        2.7, 2.7,
        3.2, 2.2,
        3.2, 1.5,
        3, 0.8,
        2.7, 0.6,
        2.5, 0,
        1.8, -0.6,
        0.5, -1.3,
        0, -1.5,
        -1.4, -1.4,
        -2, -1.2,
        -2.5, -0.8,
        -3, -0.4,
        -4, -0.5,
        -5, -0.7,
        -5.5, -0.5,
        -6, -0.3,
        -6.1, 0.3,
        -5.6, 0.2,
        -5.4, 0,
        -5, 0.3,
        -4.5, 0,
        -3, 0.5,
        -2.5, 0.4,
        -1.8, 0.5,
        -2, 1.5,
        -1.6, 1.8,
        -0.7, 2.5,
        0.3, 2.5;

    C.col(0).array() += x;
    C.col(1).array() += y;
    
    gsBSpline<T> *B =new gsBSpline<T>(kv, give(C));
    B->reverse();
    return B;
}


template<class T> gsBSpline<T> *
gsNurbsCreator<T>::BSplineFish(T const & r, T const & x, T const & y)
{
    gsKnotVector<T> kv(0,1,13,3,1,2);
    gsMatrix<T> C(16,2);

    C << -3,0,
        0,-3,
        2,-3,
        2.5,-2.5,
        4,-2.2,
        6,-2,
        6.5,-1,
        5.5,-0.5,
        5.5, 0.5,
        6,1,
        5.5,2,
        4.2,2,
        3,1.5,
        2,3,
        0,3,
        -3,0;
    C.col(0).array() += x;
    C.col(1).array() += y;
    
    gsBSpline<T> *B =new gsBSpline<T>(kv, give(C));
    return B;
}

template<class T> gsBSpline<T> *
gsNurbsCreator<T>::BSplineAmoeba3degree(T const & r, T const & x, T const & y)
{
    gsKnotVector<T> kv(0,1,18,4,1,3);
    gsMatrix<T> C(22,2);


    C << -5, -1,
        -5,0.5,
        -4.5, 1,
        -3.4, 1.05,
        -3.5, 2,
        -2.5, 3.25,
        -1.5, 3,
        -0.25, 1.5,
        0.25 ,1.6,
        2, 3,
        3.5, 3,
        4, 1.5,
        2.5, 0,
        1.9 ,-0.5 ,
        2.5, -1,
        3, -3,
        1, -4,
        -1, -3,
        -1, -1,
        -3, -1.5,
        -5, -1.5,
        -5, -1;

    C.col(0).array() += x;
    C.col(1).array() += y;
    
    gsBSpline<T> *B =new gsBSpline<T>(kv, give(C));
    B->reverse();
    return B;
}


template<class T> gsTensorNurbs<2,T> *
gsNurbsCreator<T>::NurbsDisk(T const & r, T const & x, T const & y)
{
    gsKnotVector<T> kv(0,1,0,3);
    gsMatrix<T> C(9,2);

    C << 0, -2 ,  2,-2 , 2, 0
        ,-2, -2 ,  0, 0 , 2 ,2
        ,-2 , 0 , -2, 2 , 0 ,2 ;

    C *= r;

    C.col(0).array() += x;
    C.col(1).array() += y;
    
    gsMatrix<T> ww(9, 1 ) ;
    ww(0)= 1;
    ww(1)= 0.707106781186548 ;
    ww(2)= 1;
    ww(3)= 0.707106781186548 ;
    ww(4)= 0.5 ;
    ww(5)= 0.707106781186548 ;
    ww(6)= 1;
    ww(7)= 0.707106781186548 ;
    ww(8)= 1 ;

    return new gsTensorNurbs<2,T>(kv,kv, give(C), give(ww));
}


template<class T> gsTensorBSpline<2,T> * 
gsNurbsCreator<T>::NurbsQrtPlateWHoleC0()
{
    // SK:
    // NOTE: This is still a B-SPLINE plate-with-hole!!!
    // The circular part is NOT exact.
    // TODO: Make a real NURBS geometry out of this.

    gsKnotVector<T> kv1(0,1,1,3,2);
    gsKnotVector<T> kv2(0,1,0,3);
    gsMatrix<T> C(15,2);

    C << -1, 0,
        -1, sqrt(2.0)-1,
        -1/sqrt(2.0), 1/sqrt(2.0),
        1-sqrt(2.0), 1,
        0, 1,
        -2.5, 0,
        -2.5, 0.75,
        -1.5, 1.5,
        -0.75, 2.5,
        0, 2.5,
        -4, 0,
        -4, 2,
        -4, 4,
        -2, 4,
        0, 4;

    //return new gsTensorNurbs<2,T>(kv,kv, give(C), give(ww));
    return new gsTensorBSpline<2,T>(kv1,kv2, give(C));

}


}; // namespace gismo
