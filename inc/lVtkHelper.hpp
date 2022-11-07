/**
 * @file lVtkHelper.hpp
 * @author your name (you@domain.com)
 * @brief helper functions to convert triangle data form lib3mf to vtk polydata. Also lateron prbly "pure" vtk helpers
 *  !!! This  does not use the efficient template like programming style mentioned in https://vtk.org/Wiki/VTK/Tutorials/DataArrays  !!!
 *  -> only double gets used
 * @version 0.1
 * @date 2022-03-11
 *
 * @copyright Copyright (c) 2022
 *
 */
#ifndef L_VTKHELPER_HPP
#define L_VTKHELPER_HPP

#include <vector>
#include <algorithm>
#include <iterator>

#include "lib3mf_implicit.hpp"
#include <vtkPolyData.h>
#include <vtkCell.h>
#include <vtkType.h>
#include <vtkPoints.h>
#include <vtkTriangle.h>
#include <vtkVector.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkPlaneCollection.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>

#include "VHACD.h"
#include "typesConstantsDefinition.hpp"

namespace lVtkHelper
{

    typedef double VHACD_POINT_TYPE;
    typedef uint32_t VHACD_INDEX_TYPE;

    /**
     * @brief convert 3mf mesh object to vtk poly data mesh object
     * (3mf)verts = points(vtk)
     * 
     * @param vMshObjs 
     * @return vtkSmartPointer<vtkPolyData> 
     */
    vtkSmartPointer<vtkPolyData> convertMshs2vtkPolydata(std::vector<Lib3MF::PMeshObject> vMshObjs);
    /**
     * @brief Get the Clean Polydata object / apply vtkCleanPolydata
     * 
     * @param polydata 
     * @param tolerance 
     * @return vtkSmartPointer<vtkPolyData> 
     */
    vtkSmartPointer<vtkPolyData> getCleanPolydata(vtkSmartPointer<vtkPolyData> polydata, float tolerance = 0.0);
    void calcualteCellNormalsInPlace(vtkSmartPointer<vtkPolyData> polydata, int featureAngle = 60);
    /**
     * @brief inserts the points and triangle idxs from polydata into the vectors (push_back)
     * used vectors as input instead of a struct as this allows to add mutliple data to the vectors and the advantage of a
     * struct is in this case not huge if at all (organisation wise ... for me)
     *
     * @param[in] poly  input poly data from which to get the points and tri conneciton data
     * @param[out] points output(inout) vector of points {x, y, z, x, y, z, ...} of type VHACD_POINT_TYPE (prbly. double)
     * @param[out] triangles output(inout) vector of triangle data (idxs) {idxV0, idxV1, idxV2, ...} of type VHACD_INDEX_TYPE (prbly. uint32_t)
     */
    void convertVtkPolyDataToFlatVectors(vtkPolyData *poly,
                                         std::vector<VHACD_POINT_TYPE> &vPoints, std::vector<VHACD_INDEX_TYPE> &vTriangles);

    /**
     * @brief inserts the points and triangle idxs from Lib3MF::PMeshObject into the vectors (push_back)
     *  This method should be in another file as this files name is missleading
     * @param[in] vMshObjs input vector of PMeshObjects
     * @param[out] vPoints output(inout) vector of points {x, y, z, x, y, z, ...} of type VHACD_POINT_TYPE (prbly. double)
     * @param[out] vTriangles output(inout) vector of triangle data (idxs) {idxV0, idxV1, idxV2, ...} of type VHACD_INDEX_TYPE (prbly. uint32_t)
     */
    void convertL3MFVecMeshObjsToFlatVectors(std::vector<Lib3MF::PMeshObject> vMshObjs,
                                             std::vector<VHACD_POINT_TYPE> &vPoints, std::vector<VHACD_INDEX_TYPE> &vTriangles);

    /**
     * @brief converts a vhacd convex hull to poly data. If the hull is empty successfull gets assigned false, else true
     * 
     * @param ch 
     * @param successfull 
     * @return vtkSmartPointer<vtkPolyData> 
     */
    vtkSmartPointer<vtkPolyData> convertVHACDCH2Polydata(VHACD::IVHACD::ConvexHull ch, bool &successfull);
    /**
     * @brief converts a vhacd convex hull to poly data. If the hull is empty successfull gets assigned false, else true
     * 
     * @param ch 
     * @return vtkSmartPointer<vtkPolyData> 
     */
    inline vtkSmartPointer<vtkPolyData> convertVHACDCH2Polydata(VHACD::IVHACD::ConvexHull ch) {
        bool successfull;
        return convertVHACDCH2Polydata(ch, successfull);
    }

    /**
     * @brief store vtkPolyData object as obj file
     * 
     * @param obj 
     * @param filename 
     */
    void writeVtkToObj(vtkDataObject *obj, const char *filename);
    /**
     * @brief store vtkPolyData object as vtk file
     * 
     * @param obj 
     * @param filename 
     */
    void writeVtkToVtk(vtkDataObject *obj, const char *filename);

    /**
     * @brief calculates the cos(angle) for the given vectors (dot product) "angle between"
     * calculate the cos angle of normalized vectors (dot porduct) [0,1]
     *
     * @param v0 vector 0 (assumption: is normalized)
     * @param v1 vector 1 (assumption: is normalized)
     * @return double the cosine of the angle between the vectors range: [-1, 1]
     */
    inline double cosAngleNormVectors(vtkVector3d v0, vtkVector3d v1)
    {
        // https://www.omnicalculator.com/math/angle-between-two-vectors#angle-between-two-vectors-formulas
        // https://www.wikihow.com/Find-the-Angle-Between-Two-Vectors
        return v0.Dot(v1);
    }
    /**
     * @brief calcualte the angle between two normalized vectors in rad [0, pi]
     * 
     * @param v0 
     * @param v1 
     * @return double 
     */
    inline double angleNormVectors(vtkVector3d v0, vtkVector3d v1)
    {
        return acosf64(cosAngleNormVectors(v0, v1));
    }

    /**
     * @brief multiply vector by -1 (inverting the direction)
     * 
     * @tparam T 
     * @param vec 
     * @return vtkVector3<T> 
     */
    template <typename T>
    inline vtkVector3<T> invertVector(vtkVector3<T> vec) { return vtkVector3<T>(-vec.GetX(), -vec.GetY(), -vec.GetZ()); }
    // something like this exists in vtk but ?needs init of a tri??
    inline double computeTriArea(vtkPoints *pts)
    {
        double p0[3];
        double p1[3];
        double p2[3];
        pts->GetPoint(0, p0);
        pts->GetPoint(1, p1);
        pts->GetPoint(2, p2);
        return vtkTriangle::TriangleArea(p0, p1, p2);
    }

    /**
     * @brief compute the area of a triangle
     * 
     * @param cell 
     * @return double 
     */
    inline double computeTriArea(vtkCell *cell)
    {
        auto pts = cell->GetPoints();
        return computeTriArea(pts);
    }

    /**
     * @brief get the id of the greatest cell (tri) and the size
     * 
     * @param poly 
     * @param id 
     * @param size 
     */
    void getGreatestMeshCellAndSize(vtkPolyData* poly, vtkIdType& id, double& size);

    /**
     * @brief add an array containing the triangle areas to the polydata, returns the sum of all areas works inplace (adds to poly)
     * 
     * @param poly 
     * @param newArrName name of the array that gets added to the polydata
     * @return double the sum of all areas
     */
    double appendCellAreasArrayIp(vtkPolyData* poly, const char* newArrName=TRI_AREA_ARR_NAME);

    /**
     * @brief add an array containing the triangle areas to a copy of the polydata, does not work inplace, returns the data with the appended array
     * 
     * @param poly 
     * @return vtkSmartPointer<vtkPolyData> 
     */
    vtkSmartPointer<vtkPolyData> appendCellCenters(vtkPolyData* poly);

    /**
     * @brief create a polydata point data array and assign the value of the cell it belongs to (the value will be the one from the last cell that it belongs to as a point belongs to 3 tris)
     * 
     * @param poly 
     * @param cellArrayName 
     * @param pointArrayName 
     */
    void convertCellData2PointDataIp(vtkSmartPointer<vtkPolyData> poly, const char* cellArrayName, const char* pointArrayName);

    /**
     * @brief update polydata cell data from a point data array by assigning the value of the point to the tri it belongs to, the value of the last point which belongs / is checked will be kept for the cell
     * 
     * @param poly 
     * @param cellArrayName 
     * @param pointArrayName 
     */
    void updateCellDataWithPointDataIp(vtkSmartPointer<vtkPolyData> poly, const char* cellArrayName, const char* pointArrayName);

    /**
     * @brief reduce the number of cutting planes by removing cutting planes that are at least nearly the same, by some small error margin
     * 
     * @param pts 
     * @param planc 
     * @param newPlanc 
     */
    void filterPlaneCol4Cut(vtkDataArray* pts, vtkSmartPointer<vtkPlaneCollection> planc, vtkSmartPointer<vtkPlaneCollection> newPlanc);

    /**
     * @brief clip closed surface and add an array called ParentCellId to identify to which cell a cell belonged before the cut, then set labels to 1 if they were 1 before the cut (a cutting interface)
     * 
     * @param poly 
     * @param planeCol 
     * @return vtkSmartPointer<vtkPolyData> 
     */
    vtkSmartPointer<vtkPolyData> clipClosedKeepCellLabels(vtkSmartPointer<vtkPolyData> poly, vtkSmartPointer<vtkPlaneCollection> planeCol);

    /**
     * @brief return if a rotation (by 90°) will align the printer bbx's longest side (x, y) with the partBbx's longest side (x, y)
     * 
     * @param printerBbx 
     * @param partBbx 
     * @return double the angle (0°, 90°)
     */
    double getAlignBbxsAngle(vtkBoundingBox printerBbx, vtkBoundingBox partBbx);

    /**
     * @brief create a 4x4 csv matrix from the array of size 16 representing a 4x4 matrix. ",": number separation / column, ";\n": line separation / row
     * 
     * @param vals 
     * @return std::string 
     */
    std::string getStringFromvtk4x4MatrixValues(double vals[16]);

    /**
     * @brief create a 4x4 csv matrix from the vtkTransform representing a 4x4 matrix. ",": number separation / column, ";\n": line separation / row
     * 
     * @param transf 
     * @return std::string 
     */
    inline std::string getStringFromvtk4x4MatrixValues(vtkSmartPointer<vtkTransform> transf) {
        double matVals[16];
        // copy array values, deepcopy told id would be deprec
        auto a = transf->GetMatrix()->GetData();
        std::copy(a, a+15, std::begin(matVals));
        return getStringFromvtk4x4MatrixValues(matVals);
    }

    /**
     * @brief print a 4x4 csv matrix from the vtkTransform representing a 4x4 matrix. ",": number separation / column, ";\n": line separation / row
     * 
     * @param transf 
     * @return std::string 
     */
    inline std::string printStringFromvtk4x4MatrixValues(vtkSmartPointer<vtkTransform> transf) {
        std::string theStr = getStringFromvtk4x4MatrixValues(transf);
        std::cout << theStr << std::endl;
        return theStr;
    }

}

#endif // L_VTKHELPER_HPP