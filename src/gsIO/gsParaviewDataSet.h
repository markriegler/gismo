/** @file gsParaviewDataSet.h

    @brief Provides a helper class to write Paraview (.vts) files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Karampatzakis, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsMSplines/gsMappedBasis.h>   // Only to make linker happy
#include <gsCore/gsDofMapper.h>         // Only to make linker happy
#include <gsAssembler/gsExprHelper.h>  
#include <gsAssembler/gsExprEvaluator.h>

#include<fstream>

namespace gismo 
{

class GISMO_EXPORT gsParaviewDataSet // a collection of .vts files 
{
private:
    std::string m_basename;
    gsExprHelper<real_t>::geometryMap * m_geoMap;
    gsExprEvaluator<real_t> * m_evaltr;
    index_t m_numPatches;
    gsOptionList m_options;
    
public:
    gsParaviewDataSet(std::string basename, 
                      gsExprHelper<real_t>::geometryMap * geoMap,
                      gsExprEvaluator<real_t> * eval, gsOptionList options);

                       


    template<class E>
    void addField(const expr::_expr<E> & expr,
                  std::string label)
    {
        // evaluates the expression and appends it to the vts files
        //for every patch
        unsigned nPts = m_options.askInt("numPoints",1000);
        unsigned precision = m_options.askInt("precision",5);

        std::vector<std::string> tags = m_evaltr->expr2vtk(expr, label,nPts,precision);
        std::vector<std::string> fnames = filenames();

        for ( index_t k=0; k!=m_numPatches; k++) // For every patch.
        {
            std::ofstream file;
            file.open( fnames[k].c_str(), std::ios_base::app); // Append to file
            file << tags[k];
            file.close(); 
        }
    }

    // Just here to stop the recursion
    void addFields(std::vector<std::string> labels){} 


    // The recursive case: we take a number, alongside
    // some other numbers, and produce their sum.
    template <class E, typename... Rest>
    void addFields(std::vector<std::string> labels, const expr::_expr<E> & expr, Rest... rest) {
        // keep all but first label 
        std::vector<std::string> newlabels(labels.cbegin()+1, labels.cend());
        

        addField(   expr, labels[0]);       // Add the expression 'expr' with it's corresponding label ( first one )
        addFields(   newlabels, rest...);   // Recursion
    }

    std::vector<std::string> filenames();

    void save();

    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addInt("numPoints", "Number of points per-patch.", 1000);
        opt.addInt("precision", "Number of decimal digits.", 5);
        return opt;
    }

    gsOptionList & options() {return m_options;}

};
} // End namespace gismo
