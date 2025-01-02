/** @file scaleTo01.cpp

    @brief Pre-processing step for gsHLBFGS experiments: scale the
    data so that they fit to [0, 1]^3.

	TODO: Make the inverted option available as well.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokriš
*/

#include <gismo.h>
#include <gsModeling/gsModelingUtils.hpp>

using namespace gismo;

int main(int argc, char *argv[])
{
	// Parsing the cmd arguments.
	std::string fin("");
	std::string fout("");

	index_t uvIdIn = 0;
	index_t uvIdOut = 1;
	index_t xyzIdIn = 1;
	index_t xyzIdOut = 0;

	// Intentionally negative interval to see if options are set.
	real_t tMin = 1;
	real_t tMax = 0;

	bool geo = false;
	bool verbose = false;

	gsCmdLine cmd("Scaling data to/from [0, 1]^3.\nWhen setting m and n, the data are scaled FROM [0, 1]^d, otherwise they are scaled TO [0, 1]^d.");
	cmd.addString("i", "fin", "filename of the input", fin);
	cmd.addString("o", "fout", "filename of the output", fout);

	cmd.addInt("u", "uvIdIn",   "id of the input matrix with the uv coordinates",   uvIdIn);
	cmd.addInt("v", "uvIdOut",  "id of the output matrix with the uv coordinates",  uvIdOut);
	cmd.addInt("x", "xyzIdIn",  "id of the input matrix with the xyz coordinates",  xyzIdIn);
	cmd.addInt("y", "xyzIdOut", "id of the output matrix with the xyz coordinates", xyzIdOut);

	cmd.addReal("m", "tMin", "required minimum coordinate", tMin);
	cmd.addReal("n", "tMax", "required maximum coordinate", tMax);

	cmd.addSwitch("g", "geometry", "if set to true, a geometry is scaled, otherwise a point cloud", geo);
	cmd.addSwitch("w", "verbose", "print information", verbose);

	try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

	if(geo)
		scaleGeo(fin, fout, tMin, tMax, verbose);
	else
		scalePts(fin, fout, uvIdIn, uvIdOut, xyzIdIn, xyzIdOut, tMin, tMax, verbose);

	return 0;
}
