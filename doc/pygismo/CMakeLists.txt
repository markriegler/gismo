######################################################################
## CMakeLists.txt --- doc
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

project(pyGismoDoc)

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

#find_package(sphinx)
#find_package(sphinx_rtd_theme)

# check pygismo enabled

set(SPHINX_FOUND TRUE)

if(SPHINX_FOUND)

  add_custom_target(
    doc-pygismo
    COMMAND sphinx-build -b html -d ${CMAKE_CURRENT_BINARY_DIR}/doctrees . ${CMAKE_CURRENT_BINARY_DIR}/html
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating pygismo-doc: ${gismo_BINARY_DIR}/doc/pygismo/html/index.html" VERBATIM
    )

  add_dependencies(doc doc-pygismo)

else (SPHINX_FOUND)

  add_custom_target(doc-pygismo
    ${CMAKE_COMMAND} -E cmake_echo_color --cyan "Documentation for pygismo cannot be compiled! Please install sphinx, sphinx_rtd_theme -- http://sphinx-doc.org , or python3-sphinx in package manager, or sphinx in pip -- and enable GISMO_WITH_PYBIND11 to build the pygismo  documentation."
#    ${CMAKE_COMMAND} -E cmake_echo_color --cyan "Documentation for pygismo cannot be built"
    )

endif(SPHINX_FOUND)
