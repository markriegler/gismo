
#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsBSplineBasis.h>

#include <gsNurbs/gsCompactKnotVector.h>

#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsTensorBSplineBasis.hpp>

namespace gismo
{

// instantiate tensor B-spline bases
//CLASS_TEMPLATE_INST gsTensorBasis1D<   gsBSplineBasis<real_t, gsKnotVector<real_t> > >;
CLASS_TEMPLATE_INST gsTensorBSplineBasis<1, real_t, gsKnotVector<real_t> >;
CLASS_TEMPLATE_INST gsTensorBSplineBasis<2, real_t, gsKnotVector<real_t> >;
CLASS_TEMPLATE_INST gsTensorBSplineBasis<3, real_t, gsKnotVector<real_t> >;
CLASS_TEMPLATE_INST gsTensorBSplineBasis<4, real_t, gsKnotVector<real_t> >;

// instantiate tensor B-spline bases with compact knot vector
//CLASS_TEMPLATE_INST gsTensorBasis1D< gsBSplineBasis<real_t, gsCompactKnotVector<real_t> > >;
CLASS_TEMPLATE_INST gsTensorBSplineBasis<1, real_t, gsCompactKnotVector<real_t> >;
CLASS_TEMPLATE_INST gsTensorBSplineBasis<2, real_t, gsCompactKnotVector<real_t> >;
CLASS_TEMPLATE_INST gsTensorBSplineBasis<3, real_t, gsCompactKnotVector<real_t> >;
CLASS_TEMPLATE_INST gsTensorBSplineBasis<4, real_t, gsCompactKnotVector<real_t> >;


CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<1,real_t, gsKnotVector<real_t> > >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<2,real_t, gsKnotVector<real_t> > >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<3,real_t, gsKnotVector<real_t> > >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<4,real_t, gsKnotVector<real_t> > >;

CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<1,real_t, gsCompactKnotVector<real_t> > >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<2,real_t, gsCompactKnotVector<real_t> > >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<3,real_t, gsCompactKnotVector<real_t> > >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<4,real_t, gsCompactKnotVector<real_t> > >;

}
