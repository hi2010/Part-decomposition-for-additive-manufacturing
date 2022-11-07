#ifndef L_FILEHELPER_HPP

#include "l3mfHelper.hpp"
#include "lVtkHelper.hpp"

#include <lib3mf_implicit.hpp>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkFillHolesFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkNew.h>
#include <vtkPolyDataConnectivityFilter.h>

#include <string>

class CFileHelper {
    public:
    std::string lastFp;
    Lib3MF::PModel lastPModel;

    /**
     * @brief load a 3mf mesh file to a vtkPolyData object
     * 
     * @param filePath 
     * @return vtkSmartPointer<vtkPolyData> 
     */
    vtkSmartPointer<vtkPolyData> load3mfMeshAsVtk(std::string filePath) {
        auto pmod = l3mfHelper::getPModelFrom3mfFile(filePath);
        std::vector<Lib3MF::PMeshObject> meshObjs = l3mfHelper::getMeshObjsFrom3mfModel(pmod);
        auto newMesh = vtkSmartPointer<vtkPolyData>::New();
        newMesh = lVtkHelper::convertMshs2vtkPolydata(meshObjs);
        this->lastFp = filePath;
        this->lastPModel = pmod;
        return newMesh;
    }

    /**
     * @brief try to clean, repair triangle poly data and fill holes
     * 
     * @param poly 
     * @return vtkSmartPointer<vtkPolyData> the cleaned and closed poly data if successfull
     */
    vtkSmartPointer<vtkPolyData> tryCleanAndRepairPolyData(vtkSmartPointer<vtkPolyData> poly) {
        vtkNew<vtkFillHolesFilter> fillHolesFtlr;
        vtkNew<vtkGeometryFilter> geomFltr;
        fillHolesFtlr->SetHoleSize(1000);
        auto intermPoly = lVtkHelper::getCleanPolydata(poly);
        fillHolesFtlr->SetInputData(intermPoly);
        fillHolesFtlr->Update();
        geomFltr->SetInputData(lVtkHelper::getCleanPolydata(fillHolesFtlr->GetOutput()));
        geomFltr->Update();
        auto surfPoly = geomFltr->GetOutput();
        auto cleanPoly = vtkSmartPointer<vtkPolyData>::New();
        cleanPoly = lVtkHelper::getCleanPolydata(surfPoly);
        lVtkHelper::calcualteCellNormalsInPlace(cleanPoly);
        return cleanPoly;
    }

    /**
     * @brief check how many defects get detected, NonManifold, BoundaryEdges, unconnected regions, ...
     * 
     * @param poly 
     * @return int number of defect cells and edges
     */
    int checkForDefects(vtkSmartPointer<vtkPolyData> poly) {
        // check for invalid edges
        vtkNew<vtkFeatureEdges> featureEdges;
        featureEdges->FeatureEdgesOff();
        featureEdges->BoundaryEdgesOn();
        featureEdges->NonManifoldEdgesOn();
        featureEdges->SetInputData(poly);
        int nDefectEdges = featureEdges->GetOutput()->GetNumberOfCells();

        // check if there are some connectivity errors
        // idk r8now what this is but it occured, maybe some missing indices?
        vtkNew<vtkPolyDataConnectivityFilter> connFilter;
        connFilter->SetInputData(poly);
        connFilter->SetExtractionModeToAllRegions();
        connFilter->Update();
        auto nRegions = connFilter->GetNumberOfExtractedRegions() - 1;

        auto totalDefects = nDefectEdges + nRegions;

        return totalDefects;
    }

    /**
     * @brief load 3mf mesh file to a vtkPolyData object and try to fix mesh defects. Assumption: manifold closed structure
     * 
     * @param filePath 
     * @param forceRepair try repairing even if no errors got detected
     * @return vtkSmartPointer<vtkPolyData> 
     */
    vtkSmartPointer<vtkPolyData> load3mfMeshAsVtkAndTryRepair(std::string filePath, bool forceRepair=false) {
        auto newMesh = vtkSmartPointer<vtkPolyData>::New();
        newMesh = this->load3mfMeshAsVtk(filePath);
        auto nDefects = this->checkForDefects(newMesh);
        std::cout << "input file had " << nDefects << " defects (unconnected regions, and defect edges summed up)" << std::endl;
        if (nDefects == 0 && forceRepair == false) {
            return newMesh;
        }
        std::cout << "trying to repair input file" << std::endl;
        newMesh = this->tryCleanAndRepairPolyData(newMesh);
        nDefects = this->checkForDefects(newMesh);
        std::cout << "after repair had " << nDefects << " defects (unconnected regions, and defect edges summed up)" << std::endl;
        return newMesh;
    }

};
#define L_FILEHELPER_HPP
#endif //L_FILEHELPER_HPP