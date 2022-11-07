#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>
#include <vtkCleanPolyData.h>
#include <vtkCellIterator.h>
#include <vtkOBJWriter.h>
#include <vtkAppendLocationAttributes.h>
#include <vtkCell.h>
#include <vtkPolyDataWriter.h>
#include <vtkStaticCleanPolyData.h>
#include <vtkClipClosedSurface.h>
#include <vtkAssume.h>
#include <vtkArrayDispatch.h>
#include <vtkDataArrayAccessor.h>
#include <vtkType.h>
#include <vtkAssume.h>

#include "lVtkHelper.hpp"
#include "arrayMath.hpp"
#include "modifiedClipClosedSurface.hpp"

namespace lVtkHelper
{
    // if true: some functions use some asserts that are needed for this project
    static const bool USE_ASSERTS = true;
    typedef vtkDoubleArray LVH_POINTS_DATATYPE;
    // gets automatically assigned float or double, depending on LVH_POINTS_DATATYPE
    typedef std::conditional<std::is_same<LVH_POINTS_DATATYPE, vtkDoubleArray>::value, double, float>::type LVH_POINT_DATATYPE;

    vtkSmartPointer<vtkPolyData> convertMshs2vtkPolydata(std::vector<Lib3MF::PMeshObject> vMshObjs)
    {
        auto polydata = vtkSmartPointer<vtkPolyData>::New();
        vtkNew<vtkPoints> points;
        vtkNew<vtkCellArray> polys;
        // load 3mf model into vtkPolydatatype and casting data to double
        for (auto &obj : vMshObjs)
        {
            // add the points to the points data
            std::vector<Lib3MF::sPosition> verts;
            obj->GetVertices(verts);
            for (auto &v : verts)
            {
                // cast to double to get more precise calculations in the vtk calculations -> needs to be back converted when new parts get created
                LVH_POINT_DATATYPE p[3] = {static_cast<LVH_POINT_DATATYPE>(v.m_Coordinates[0]), static_cast<LVH_POINT_DATATYPE>(v.m_Coordinates[1]), static_cast<LVH_POINT_DATATYPE>(v.m_Coordinates[2])};
                points->InsertNextPoint(p);
            }
            //auto triCount = obj->GetTriangleCount();
            std::vector<Lib3MF::sTriangle> tris;
            obj->GetTriangleIndices(tris);
            // add the tri idxs to the polys data
            for (auto &tri : tris)
            {
                polys->InsertNextCell(3);
                // polys->InsertNextCell(VTK_TRIANGLE);
                polys->InsertCellPoint(tri.m_Indices[0]);
                polys->InsertCellPoint(tri.m_Indices[1]);
                polys->InsertCellPoint(tri.m_Indices[2]);
            }
        }
        // assing points to the polydata
        polydata->SetPoints(points);
        // assing polys / cells / the connection information to the polydata
        polydata->SetPolys(polys);

        return polydata;
    }

    vtkSmartPointer<vtkPolyData> getCleanPolydata(vtkSmartPointer<vtkPolyData> polydata, float tolerance /*= 0*/)
    {
        // fixup poly data (remove duplicates and change degenrates)
        vtkNew<vtkStaticCleanPolyData> cleanPolyData;
        // tolerance == 0.0 changes the used algo and is faster
        cleanPolyData->SetTolerance(tolerance);
        cleanPolyData->SetInputData(polydata);
        cleanPolyData->Update();
        vtkSmartPointer<vtkPolyData> cleanedPoly = vtkSmartPointer<vtkPolyData>::New();
        cleanedPoly = cleanPolyData->GetOutput();
        return cleanedPoly;
    }

    void calcualteCellNormalsInPlace(vtkSmartPointer<vtkPolyData> polydata, int featureAngle /*= 60*/)
    {
        vtkNew<vtkPolyDataNormals> surfNormals;
        surfNormals->SetInputData(polydata);
        // this is the angle for sharp corners, as on sharp corners some special splitting calc is done
        surfNormals->SetFeatureAngle(featureAngle);
        surfNormals->SetComputeCellNormals(true);
        surfNormals->SetComputePointNormals(false);
        surfNormals->Update();
        auto snOut = surfNormals->GetOutput();
        auto norms = vtkArrayDownCast<vtkFloatArray>(snOut->GetCellData()->GetArray("Normals"));
        polydata->GetCellData()->SetNormals(norms);

        if (USE_ASSERTS) {
            for (vtkIdType i = 0; i < polydata->GetNumberOfCells(); i++){
                assert(polydata->GetCellType(i) == VTK_TRIANGLE);
            }
            // crashes on empty polydata
            assert(polydata->GetNumberOfCells() == polydata->GetCellData()->GetNormals()->GetNumberOfTuples());
        }
    }

    void convertVtkPolyDataToFlatVectors(vtkPolyData *poly,
                                         std::vector<VHACD_POINT_TYPE> &vPoints, std::vector<VHACD_INDEX_TYPE> &vTriangles)
    {
        auto nPoints = poly->GetNumberOfPoints();
        auto pts = poly->GetPoints();
        for (vtkIdType i = 0; i < nPoints; i++)
        {
            auto p = pts->GetPoint(i);
            vPoints.push_back(static_cast<VHACD_POINT_TYPE>(p[0]));
            vPoints.push_back(static_cast<VHACD_POINT_TYPE>(p[1]));
            vPoints.push_back(static_cast<VHACD_POINT_TYPE>(p[2]));
        }

        vtkSmartPointer<vtkCellIterator> it = poly->NewCellIterator();
        vtkNew<vtkGenericCell> cell;
        for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell())
        {
            it->GetCell(cell);
            auto nOfPts = cell->GetNumberOfPoints();
            for (vtkIdType i = 0; i < nOfPts; i++)
            {
                vTriangles.push_back(static_cast<VHACD_INDEX_TYPE>(cell->GetPointId(i)));
            }
        }
        it->Delete();
    }

    void convertL3MFVecMeshObjsToFlatVectors(std::vector<Lib3MF::PMeshObject> vMshObjs,
                                             std::vector<VHACD_POINT_TYPE> &vPoints, std::vector<VHACD_INDEX_TYPE> &vTriangles)
    {
        // insert all mesh data in the vectors for vhacd to process
        for (auto &obj : vMshObjs)
        {
            std::vector<Lib3MF::sPosition> verts;
            obj->GetVertices(verts);
            for (auto &v : verts)
            {
                // cast to double to get more precise calculations in the vtk calculations -> needs to be back converted when new parts get created
                vPoints.push_back(static_cast<VHACD_POINT_TYPE>(v.m_Coordinates[0]));
                vPoints.push_back(static_cast<VHACD_POINT_TYPE>(v.m_Coordinates[1]));
                vPoints.push_back(static_cast<VHACD_POINT_TYPE>(v.m_Coordinates[2]));
            }
            std::vector<Lib3MF::sTriangle> tris;
            obj->GetTriangleIndices(tris);
            for (auto &tri : tris)
            {
                vTriangles.push_back(static_cast<VHACD_INDEX_TYPE>(tri.m_Indices[0]));
                vTriangles.push_back(static_cast<VHACD_INDEX_TYPE>(tri.m_Indices[1]));
                vTriangles.push_back(static_cast<VHACD_INDEX_TYPE>(tri.m_Indices[2]));
            }
        }
    }

    vtkSmartPointer<vtkPolyData> convertVHACDCH2Polydata(VHACD::IVHACD::ConvexHull ch, bool &successfull) {
        auto newPoly = vtkSmartPointer<vtkPolyData>::New();
        vtkNew<vtkPoints> pts;
        vtkNew<vtkCellArray> tris;
        // in theory this should never happen but it did and time is ... so use some "ignore" tactic
        if (ch.m_nPoints == 0 || ch.m_nTriangles == 0) {
            successfull = false;
            std::cout << "in convert ch to poly got 0 pts or 0 tri: nTri: " << ch.m_nTriangles << " nPts: " << ch.m_nPoints << std::endl;
            return newPoly;
        }
        successfull = true;

        assert(ch.m_nPoints > 0);
        assert(ch.m_nTriangles > 0);
        // according to the test file of vhacd and the SaveOBJ, nPts, and nTris is the n of pts / tris (1/3 of entries)
        for (size_t i = 0; i < ch.m_nPoints*3; i+=3) {
            pts->InsertNextPoint(ch.m_points[i], ch.m_points[i+1], ch.m_points[i+2]);
        }
        for (size_t i = 0; i < ch.m_nTriangles*3; i+=3) {
            tris->InsertNextCell(3);
            tris->InsertCellPoint(ch.m_triangles[i]);
            tris->InsertCellPoint(ch.m_triangles[i+1]);
            tris->InsertCellPoint(ch.m_triangles[i+2]);
        }
        newPoly->SetPoints(pts);
        newPoly->SetPolys(tris);
        return newPoly;
    }

    void writeVtkToObj(vtkDataObject* obj, const char* filename) {
        vtkNew<vtkOBJWriter> writer;
        writer->SetFileName(filename);
        writer->SetInputData(obj);
        writer->Write();
    }

    void writeVtkToVtk(vtkDataObject* obj, const char* filename) {
        vtkNew<vtkPolyDataWriter> polyWrit;
        // for paraview, cause with ver 5.2 it crashes (paraview 5.7)
        polyWrit->SetFileVersion(vtkPolyDataWriter::VTKFileVersion::VTK_LEGACY_READER_VERSION_4_2);
        polyWrit->SetInputData(obj);
        polyWrit->SetFileName(filename);
        polyWrit->Write();
    }

    void getGreatestMeshCellAndSize(vtkPolyData* poly, vtkIdType& id, double& size) {
        // assumed: trimesh
        auto nCells = poly->GetNumberOfCells();
        // size of the greatest face
        double maxCellSz = 0;
        // according to the paper the printing bed normal should be used here, but for xyz the printing bed normal is just the (in this case negative) printing direction (z-Axis)
        vtkIdType greatestCellId = 0;
        double greatestCellSize = 0;
        
        for (vtkIdType i = 0; i < nCells; i++) {
            // this ones not thread safe
            vtkCell* cell = poly->GetCell(i);
            double triAr = computeTriArea(cell);
            if (triAr  > maxCellSz) {
                maxCellSz = triAr;
                greatestCellId = i;
                greatestCellSize = maxCellSz;
            }
        }
        id = greatestCellId;
        size = greatestCellSize;
    }

    double appendCellAreasArrayIp(vtkPolyData* poly, const char* newArrName) {
        // would be more robust with vtkCellSizeFilter but ...
        if (poly->GetCellData()->HasArray(newArrName)) {
            poly->GetCellData()->RemoveArray(newArrName);
        }
        vtkSmartPointer<vtkFloatArray> triAreas = vtkSmartPointer<vtkFloatArray>::New();
        triAreas->SetNumberOfComponents(1);
        triAreas->SetName(newArrName);
        double totalArea = 0;
        double p0[3];
        double p1[3];
        double p2[3];
        auto nCells = poly->GetNumberOfCells();
        for (vtkIdType i=0; i < nCells; i++ )
        {
            vtkCell* tri = poly->GetCell(i);
            tri->GetPoints()->GetPoint(0, p0);
            tri->GetPoints()->GetPoint(1, p1);
            tri->GetPoints()->GetPoint(2, p2);
            // TODO: make this more safe
            // for all points are the same a high value gets returend
            double area = vtkTriangle::TriangleArea(p0, p1, p2);
            //if (abs(area) > .001) std::cout << "big ar" << area << std::endl;
            totalArea += area;
            triAreas->InsertNextValue(area);
        }
        poly->GetCellData()->AddArray(triAreas);
        return totalArea;
    }

    vtkSmartPointer<vtkPolyData> appendCellCenters(vtkPolyData* poly) {
        vtkNew<vtkAppendLocationAttributes> fltrCellCenters;
        fltrCellCenters->AppendCellCentersOn();
        fltrCellCenters->SetInputData(poly);
        fltrCellCenters->Update();
        vtkSmartPointer<vtkPolyData> res = fltrCellCenters->GetPolyDataOutput();
        return res;
    }

    void convertCellData2PointDataIp(vtkSmartPointer<vtkPolyData> poly, const char* cellArrayName, const char* pointArrayName) {
        // kiinda hack -> set for every point on an interface the point data to interface afterwards, if all cell pts are marked as if -> the cell is an if  -> enables keeping cell data
        if (poly->GetCellData()->HasArray(cellArrayName) == false) {
            return;
        }
        vtkNew<vtkFloatArray> idar;
        idar->SetName(pointArrayName);
        idar->SetNumberOfComponents(1);
        idar->SetNumberOfTuples(poly->GetNumberOfPoints());
        for (vtkIdType i = 0; i < poly->GetNumberOfCells(); i++) {
            auto cell = poly->GetCell(i);
            for (vtkIdType j = 0; j  < cell->GetNumberOfPoints(); j++){
                auto ptId = cell->GetPointId(j);
                if (poly->GetCellData()->GetArray(cellArrayName)->GetComponent(i, 0) > 0) {
                    idar->SetComponent(ptId, 0, 1);
                } else {
                    idar->SetComponent(ptId, 0, 0);
                }
            }
        }
        poly->GetPointData()->AddArray(idar);
    }

    void updateCellDataWithPointDataIp(vtkSmartPointer<vtkPolyData> poly, const char* cellArrayName, const char* pointArrayName) {
        // has been created to overwrite the labels array.
        // checks if all points of a cell have a value of ~1 or up. If so the cellArray gets updated with a 1 at that cell, else not
        if (poly->GetPointData()->HasArray(pointArrayName) == false) {
            return;
        }
        // if used and not set to simply 1, every new cut should get new colors
        // IF the interfacepoint array gets deleted after every cut -> no the array get overwritten all the time which migt be ok or not.
        //auto newValue = poly->GetCellData()->GetArray(cellArrayName)->GetMaxDiscreteValues() + 1;
        auto newValue = 1;
        for (vtkIdType i = 0; i < poly->GetNumberOfCells(); i++) {
            auto cell = poly->GetCell(i);
            auto isIf = true;
            if (cell->GetNumberOfPoints() == 0) isIf = false;
            for (vtkIdType j = 0; j  < cell->GetNumberOfPoints(); j++){
                auto ptId = cell->GetPointId(j);
                if (poly->GetPointData()->GetArray(pointArrayName)->GetComponent(ptId, 0) < .9999) {
                    isIf = false;
                    break;
                }
            }
            if (isIf) {
                poly->GetCellData()->GetArray(cellArrayName)->SetComponent(i, 0, newValue);
            }
        }
    }

    // functor for faster filtering of planes that do not reallly cut the object
    struct filterPlaneColNoCut {
        template <typename PointsArray>
        void operator()(PointsArray *pts, vtkSmartPointer<vtkPlaneCollection> planc, vtkSmartPointer<vtkPlaneCollection> newPlanc) {
            VTK_ASSUME(pts->GetNumberOfComponents() == 3);

            vtkDataArrayAccessor<PointsArray> p(pts);
            vtkIdType nPts = pts->GetNumberOfTuples();
            auto nPlanes = planc->GetNumberOfItems();
            
            for (vtkIdType j = 0; j < nPlanes; j++) {
                bool pointOutside = false;
                auto plane = planc->GetItem(j);
                for (vtkIdType i = 0; i < nPts; ++i) {
                    double pt[3] = {p.Get(i, 0), p.Get(i, 1), p.Get(i, 2)};
                    auto dist = plane->Evaluate(plane->GetNormal(), plane->GetOrigin(), pt);
                    //auto dist = 0;
                    if (dist < -std::numeric_limits<double>::min()*1000) {
                        pointOutside = true;
                        break;
                    }
                }
                if (pointOutside) {
                    newPlanc->AddItem(plane);
                }
            }
        }
    };

    void filterPlaneCol4Cut(vtkDataArray* pts, vtkSmartPointer<vtkPlaneCollection> planc, vtkSmartPointer<vtkPlaneCollection> newPlanc) {
        filterPlaneColNoCut worker;
        typedef vtkArrayDispatch::DispatchByValueType
        <
            vtkArrayDispatch::Reals
        > Dispatcher;
        if (!Dispatcher::Execute(pts, worker, planc, newPlanc)) {
            worker(pts, planc, newPlanc);
        }
    }

    vtkSmartPointer<vtkPolyData> clipClosedKeepCellLabels(vtkSmartPointer<vtkPolyData> poly, vtkSmartPointer<vtkPlaneCollection> planeCol){
        // !! Crashes if the result has 0 points !!

        if (poly->GetNumberOfCells() == 0) {
            std::cout << "gota zero cella" << std::endl;
            assert(0);
        }
        // remove planes that do not really cut anything significant (double.min() *1000) 
        // if not used there might be a chance, that cutting may cause everything to become interface area (cube)
        vtkNew<vtkPlaneCollection> newPlaneCol;
        filterPlaneCol4Cut(poly->GetPoints()->GetData(), planeCol, newPlaneCol);
        if (newPlaneCol->GetNumberOfItems() == 0) {
            return poly;
        }
        vtkNew<ModClipClosed> clipClosedFltr;

        clipClosedFltr->SetInputData(poly);
        clipClosedFltr->SetScalarModeToLabels();
        clipClosedFltr->SetClippingPlanes(newPlaneCol);
        clipClosedFltr->SetPassPointData(true);
        clipClosedFltr->GlobalWarningDisplayOn();
        clipClosedFltr->SetTriangulationErrorDisplay(true);
        clipClosedFltr->Update();
        vtkSmartPointer<vtkPolyData> newPoly = vtkSmartPointer<vtkPolyData>::New();
        newPoly = clipClosedFltr->GetOutput();
        //std::cout << "done" << std::endl;

        // update labels in the output poly, based on the input poly
        if (poly->GetCellData()->HasArray("Labels") == 0) return newPoly;

        // from here on both input and output polydata have labels
        //auto nNewCells = newPoly->GetNumberOfCells();
        auto paridAr = newPoly->GetCellData()->GetArray("ParentCellId");
        if (!paridAr) return newPoly;
        auto nRels = paridAr->GetNumberOfTuples();
        auto parLbls = poly->GetCellData()->GetArray("Labels");
        auto newLbls = newPoly->GetCellData()->GetArray("Labels");

        for (vtkIdType i = 0; i < nRels; i++){
            auto parentId = paridAr->GetComponent(i, 0);
            if (parentId == -1) continue;
            auto parVal = parLbls->GetComponent(parentId, 0);
            // label already set (everything bigger than 0.5 is considered as set (normally 1))
            if (newLbls->GetComponent(i, 0) > 0.5) continue;
            newLbls->SetComponent(i, 0, parVal);
        }

        return newPoly;
    }

    /**
     * @brief check how the longest side of the printer bbx aligns with the longes side of the partBbx / if a rotation by 90 deg is needed
     * 
     * @param printerBbx 
     * @param partBbx 
     * @return the angle to align the bbxs (0 or 90 deg)
     */
    double getAlignBbxsAngle(vtkBoundingBox printerBbx, vtkBoundingBox partBbx) {
        // check if flip might help:
        double partLens[3];
        double printLens[4];
        partBbx.GetLengths(partLens);
        printerBbx.GetLengths(printLens);
        // 0, 0 -> no rot. 1, 1 -> no rot. 1, 0 or 1, 0 rot -> if sum == 1 -> rot
        int partMaxDir = (partLens[0] > partLens[1]) ? 0 : 1;
        int printMaxDir = (printLens[0] > printLens[1]) ? 0 : 1;
        bool useRot = (partMaxDir + printMaxDir) == 1;
        // "rotate" the print bbx before cutting
        return useRot*90.;
    }

    std::string getStringFromvtk4x4MatrixValues(double vals[16]) {
        std::string assemblyStr = "";
        // to csv with sep: "," and linebreak: ";"
        for (int i = 0; i < 16; i++) {
            assemblyStr += std::to_string(vals[i]);
            if (i != 0 && (i+1) % 4 == 0) {
                // linebreak
                assemblyStr += ";\n";
            } else {
                // new value
                assemblyStr += ",";
            }
        }
        return assemblyStr;
    }

}