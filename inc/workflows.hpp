/**
 * @file workflows.hpp
 * @author your name (you@domain.com)
 * @brief the whole top level optimization workflows and data structures
 * @version 0.1
 * @date 2022-08-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef WORKFLOWS_HPP
#define WORKFLOWS_HPP

#include <string>
#include <vector>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include "VHACD.h"
#include "meshWithAttributes.hpp"
#include "openGA.hpp"

// should contain the general workflow like workflow optimize ... -> do a, b, c, .... (for one gen)

/**
 * @brief the weigths for the optimization
 * 
 */
struct sCostWeights {
    double costRoughness = 1;
    double costOverhangArea = 1;
    double costSharpness = 1;
    double costGap = 1;
    double costConcavity = 1;
    double costFeasibility = 1;
    double costInterfaceArea = 1;
    // part count
    double costQuantity = 1;

    std::vector<std::string> toStringVector(){
        std::vector<std::string> resVec;
        resVec.push_back("costRoughness"); resVec.push_back(std::to_string(costRoughness));
        resVec.push_back("costOverhangArea"); resVec.push_back(std::to_string(costOverhangArea));
        resVec.push_back("costSharpness"); resVec.push_back(std::to_string(costSharpness));
        resVec.push_back("costGap"); resVec.push_back(std::to_string(costGap));
        resVec.push_back("costConcavity"); resVec.push_back(std::to_string(costConcavity));
        resVec.push_back("costFeasibility"); resVec.push_back(std::to_string(costFeasibility));
        resVec.push_back("costInterfaceArea"); resVec.push_back(std::to_string(costInterfaceArea));
        resVec.push_back("costQuantity"); resVec.push_back(std::to_string(costQuantity));
        return resVec;
    }

    std::string toString(){
        return \
        " costRoughness: "      + std::to_string(costRoughness) + \
        " costOverhangArea: "   + std::to_string(costOverhangArea) + \
        " costSharpness: "      + std::to_string(costSharpness) + \
        " costGap: "            + std::to_string(costGap) + \
        " costConcavity: "      + std::to_string(costConcavity) + \
        " costFeasability: "    + std::to_string(costFeasibility) + \
        " costInterfaceArea: "  + std::to_string(costInterfaceArea) + \
        " costQuantity: "       + std::to_string(costQuantity);
    }
    std::string toStringCustomSep(const std::string &sepBetween="=", const std::string &sepAfter=" "){
        return \
         "costRoughness"      + sepBetween + std::to_string(costRoughness) + sepAfter + \
        " costOverhangArea"   + sepBetween + std::to_string(costOverhangArea) + sepAfter + \
        " costSharpness"      + sepBetween + std::to_string(costSharpness) + sepAfter + \
        " costGap"            + sepBetween + std::to_string(costGap) + sepAfter + \
        " costConcavity"      + sepBetween + std::to_string(costConcavity) + sepAfter + \
        " costFeasability"    + sepBetween + std::to_string(costFeasibility) + sepAfter + \
        " costInterfaceArea"  + sepBetween + std::to_string(costInterfaceArea) + sepAfter + \
        " costQuantity"       + sepBetween + std::to_string(costQuantity);
    }
    std::string toStringJson(){
        std::string resStr = "{\n";
        auto vecRep = toStringVector();
        for (size_t i = 0; i < vecRep.size()-2; i+=2) {
            resStr += '"' + vecRep.at(i) + "\": " + vecRep.at(i+1) + ",\n";
        }
        resStr += '"' + vecRep.at(vecRep.size()-2) + "\": " + vecRep.at(vecRep.size()-1) + "\n";
        resStr += std::string("}");
        return resStr;
    }
};

/**
 * @brief the costs of a decomposition
 * 
 */
struct sPartCosts: sCostWeights{
    sPartCosts(){
        costRoughness        = 0;
        costOverhangArea     = 0;
        costSharpness        = 0;
        costGap              = 0;
        costConcavity        = 0;
        costFeasibility      = 0;
        costInterfaceArea    = 0;
        costQuantity         = 0;
    }
    double totalCost = 0;

    double calculateTotalCost(){
        totalCost = 0;
        totalCost += costRoughness    ;
        totalCost += costOverhangArea ;
        totalCost += costSharpness    ;
        totalCost += costGap          ;
        totalCost += costConcavity    ;
        totalCost += costFeasibility  ;
        totalCost += costInterfaceArea;
        totalCost += costQuantity     ;
        return totalCost;
    }
};

/**
 * @brief a solution set / species for the ga
 * 
 */
struct MySolution {
     // the vhacd paramters that are varied
     double alpha, beta, gamma, maxConcavity;

     std::string to_string() const {
         return "alpha= "   + std::to_string(alpha)
         +" beta= "         + std::to_string(beta)
         +" gamma= "        + std::to_string(gamma)
         +" maxConcavity"   + std::to_string(maxConcavity);
     }
};

/**
 * @brief the middle cost of a solution
 * 
 */
struct MyMiddleCost: sPartCosts{
};

typedef EA::Genetic<MySolution, MyMiddleCost> GA_Type;
typedef EA::GenerationType<MySolution, MyMiddleCost> Generation_Type;

struct OptimizationSettings {
    struct sCostWeights costWts;
    GA_Type gaSettings;
};

/**
 * @brief calculate the costs for a single object
 * 
 * @param meshwa the object
 * @param costs the weights
 * @return double the cost
 */
double calculateCost(polyQualityMeasures::MeshWithAttributes* meshwa, sCostWeights costs);

/**
 * @brief optimize a decomposition using an ga with some default settings
 * 
 * @param meshwa the mesh for which to find the best v-hacd params
 * @return VHACD::IVHACD::Parameters the v-hacd params of the best solution
 */
VHACD::IVHACD::Parameters optimizeGA (polyQualityMeasures::MeshWithAttributes meshwa);

/**
 * @brief optimize the decomposition of a mesh
 * 
 * @param meshwa the mesh to decompose
 * @param optSettings settings for the opt algo (GA)
 * @param ptrMidCost returns the middle cost of the best result if nullptr -> ignored
 * @return VHACD::IVHACD::Parameters the best v-hacd parameters
 */
VHACD::IVHACD::Parameters optimizeGA (polyQualityMeasures::MeshWithAttributes meshwa, OptimizationSettings optSettings, MyMiddleCost* ptrMidCost=nullptr);

/**
 * @brief split object using the given params, then store the results
 *
 * @param meshwa the mesh to decompose
 * @param params the params to use
 * @param path the output path
 * @param resPrefix the resulting files prefix
 * @param hullPrefix the resulting convex hulls prefix
 * @param thisCost the cost values to store to the result folder
 */
void splitMeshwaAndSave2File(polyQualityMeasures::MeshWithAttributes* meshwa, VHACD::IVHACD::Parameters params, const std::string &path="/home/louis/CAD/debugGA/", const std::string &resPrefix="result", const std::string &hullPrefix="hullResult", MyMiddleCost* thisCost=nullptr);

#endif