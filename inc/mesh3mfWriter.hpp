#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "SimpleLogger.hpp"

// P* are pointers
#include "lib3mf_implicit.hpp"
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkNew.h>
#include <vtkTriangleFilter.h>

typedef std::vector<unsigned int>::size_type tIndIdx;

namespace mesh3mfWriter {

/**
 * @brief add the mesh data from a vtkPolyData object to a lib3mf PMeshObject
 * 
 * @param poly 
 * @param mesh3mf 
 */
inline void convertPolyData23mfMesh(vtkSmartPointer<vtkPolyData> poly, Lib3MF::PMeshObject mesh3mf) {
    // prefilter data, probly not needed but ensures that the data consists of tris
    vtkNew<vtkTriangleFilter> triFilter;
    triFilter->PassVertsOff();
    triFilter->PassLinesOff();
    triFilter->SetInputData(poly);
    triFilter->Update();
    auto triPoly = triFilter->GetOutput();

    // copy vertices / points (may contain points with no connection)
    auto nPts = triPoly->GetNumberOfPoints();
    for (vtkIdType i = 0; i < nPts; i++) {
        Lib3MF::sPosition pos;
        double pt[3];
        triPoly->GetPoint(i, pt);
        for (auto j = 0; j < 3; j++) {
            pos.m_Coordinates[j] = static_cast<Lib3MF_single>(pt[j]);
        }
        mesh3mf->AddVertex(pos);
    }

    auto nClls = triPoly->GetNumberOfCells();
    for (vtkIdType i = 0; i < nClls; i++) {
        // assumption -> all clls are tris
        vtkNew<vtkGenericCell> cll;
        triPoly->GetCell(i, cll);
        Lib3MF::sTriangle tri;
        for (vtkIdType j = 0; j < 3; j++) {
            tri.m_Indices[j] = static_cast<Lib3MF_uint32>(cll->GetPointId(j));
        }
        mesh3mf->AddTriangle(tri);
    }
    //poly->GetPoint()
}

// returns false if it failed
/**
 * @brief add a vtkPolyData objects mesh to a lib3mf PModel as new object
 * 
 * @param poly 
 * @param resModel 
 * @param wrapper 
 * @param polyName 
 * @return true always
 * @return false never
 */
inline bool addPolyData2Model(vtkSmartPointer<vtkPolyData> poly, Lib3MF::PModel resModel, Lib3MF::PWrapper wrapper, std::string polyName="objectName") {

    auto curMesh = resModel->AddMeshObject();
    convertPolyData23mfMesh(poly, curMesh);
    curMesh->SetName(polyName);
    resModel->AddBuildItem(curMesh.get(), wrapper->GetIdentityTransform());
    return true;
}

/**
 * @brief create a lib3mf PModel and add the vtkPolyData poly's mesh as new object
 * 
 * @param poly 
 * @param polyName 
 * @return Lib3MF::PModel 
 */
inline Lib3MF::PModel convertPolyData2PModel(vtkSmartPointer<vtkPolyData> poly, const std::string &polyName="objectName") {
    Lib3MF::PWrapper wrapper = Lib3MF::CWrapper::loadLibrary();
    // from here on the library loading worked
    Lib3MF::PModel model = wrapper->CreateModel();
    addPolyData2Model(poly, model, wrapper, polyName);
    return model;
}

/**
 * @brief write vtkPolyData mesh as 3mf file containing the mesh as a mesh object
 * 
 * @param poly 
 * @param outFileP 
 * @param polyName 
 * @return true alway
 * @return false never
 */
inline bool writeVtkPolyDataAs3mfModel(vtkSmartPointer<vtkPolyData> poly, std::string outFileP, const std::string &polyName="objectName") {
    
    //auto loadSuccess = loadObjFileAndAdd2Model(inFileP, model, wrapper);
    auto model = convertPolyData2PModel(poly, polyName);
    
    //logger.log(0, "3mf model created, preparing to write");
    Lib3MF::PWriter writer = model->QueryWriter("3mf");
    logger.log(0, "writer object initialized");
    writer->WriteToFile(outFileP);
    logger.log(4, "3mf file has been written to:\n" + outFileP);    
    return true;
}

}