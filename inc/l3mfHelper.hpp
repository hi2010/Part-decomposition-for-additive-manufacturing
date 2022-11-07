/**
 * @file l3mfHelper.hpp
 * @author your name (you@domain.com)
 * @brief helper for reading of 3mf files
 * @version 0.1
 * @date 2022-08-24
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef L3MF_HELPER_H
#define L3MF_HELPER_H

#include <string>

#include "lib3mf_implicit.hpp"

namespace l3mfHelper{

    /**
     * @brief print warnings from the reader (after reading a file)
     * 
     * @param reader
     */
    void printReaderWarnings(Lib3MF::PReader reader);
    Lib3MF::PModel getPModelFrom3mfFile(std::string filename, bool bPrintReaderWarnings=true);
    void printMetadataGroup(Lib3MF::PMetaDataGroup meta);
    /**
     * @brief print the version of lib3mf
     * 
     * @param wrapper 
     */
    void print3mfVersion(Lib3MF::PWrapper wrapper);
    /**
     * @brief Get all mesh objects from the 3mf model object
     * 
     * @param model 
     * @return std::vector<Lib3MF::PMeshObject> 
     */
    std::vector<Lib3MF::PMeshObject> getMeshObjsFrom3mfModel(Lib3MF::PModel model);

}
#endif // L3MF_HELPER_H