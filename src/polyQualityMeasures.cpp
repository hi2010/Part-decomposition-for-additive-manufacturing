//#include <stdexcept>

#include <math.h>
#include <limits>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkTriangle.h>
#include <vtkMath.h>
#include <vtkTuple.h>
#include <vtkMassProperties.h>
#include <vtkFloatArray.h>
#include <vtkAbstractCellLocator.h>
#include <vtkModifiedBSPTree.h>
#include <vtkMapper.h>
#include <vtkTriangleFilter.h>
#include <vtkMultiObjectMassProperties.h>
#include <vtkLine.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkMath.h>

#include "arrayMath.hpp"
#include "polyQualityMeasures.hpp"
#include "typesConstantsDefinition.hpp"

namespace polyQualityMeasures {
    // TODO: move the comments to the header

    /* #region quality measures as described in Oh et al. Part Decomposition and Evaluation ... */
    //---------------------------------------------------------------------------------------------------

    /**
     * @brief Calculates the cost funciton for the roughness of one component
     * 
     * @param meshwa -> mesh: the polydata for which to calculate 
     *  assumption: normal vectors are calculated and normalized otherwise i would need to check it every time which is kinda not needed in this case
     * ->printDir the build direction of the printer (will get normalized)
     * @return double the calculated cost
     */
    double roughness(MeshWithAttributes* meshwa) {
        // the 4 here is kinda rand. If on assumes approx .5 min error per calculation this means upto 8 erroneous calcuations
        // but this should be taken with a grain of salt as there could be worse cases so... jus a rndm assumption
        auto errRange = std::numeric_limits<double>::min()*4;
        auto upperErrRange = 1 - errRange;
        auto buildDir = meshwa->printDir;
        buildDir.Normalize();

        // iterate over the faces of the component to calculate the per face cost and add it up as component cost
        double totalCost = 0.0;
        auto cellNorms = meshwa->getMeshCellNormals();
        auto cellAreas = meshwa->getMeshTriAreas();
        auto nNorms = cellNorms->GetNumberOfTuples();
        for (vtkIdType i = 0; i < nNorms; i++) {
            auto vNorm = vtkVector3d(cellNorms->GetTuple(i));
            // vNorm.normalize(); should not be used as i am calculating normalized normals using vtkPolyDataNormals for calculation
            // the cost function from the paper says the angle should lie in between (0, pi). For cosAngle this is (-1,1), otherwise use 0
            // this is stange cause the angle is inbetween [0, pi] which would make the value at 0 and pi == 0 and jump in the direct
            // vincinity to 1. On the other hand this seems legit, cause for angle = 0 and pi this means the face is horizontal
            // -> any small deviation means the face is very flat -> hich staristep effect
            // for reasons of calculation errors, ... the absolute 0 and pi should be converted to a small range
            // this converts to: if cosAngle in [-e,e] or in [1-e, 1+e] := 0
            auto diffAngle = lVtkHelper::cosAngleNormVectors(vNorm, buildDir);
            // !!! there is a not at the begining
            // err range is because rounding error might avoid landingn on 0 or pi
            if (!((diffAngle > errRange) || (diffAngle < upperErrRange))) {
                // the angle is within the range
                auto absDiffAngle = fabs(diffAngle);
                double triAr = cellAreas->GetComponent(i,0);
                totalCost += absDiffAngle * triAr;
            }
        }

        if (printResults) std::cout << __func__ << " run with points: " << totalCost << std::endl;
        return totalCost;
    }

    double overhangArea(MeshWithAttributes* meshwa) {
        // assumed: trimesh
        auto poly = meshwa->mesh;
        auto nCells = poly->GetNumberOfCells();
        auto cellAreas = meshwa->getMeshTriAreas();
        // critical overhang angle in rad
        double overhangAngle = meshwa->overhangAngle;
        double criticalCosAngle = cos(overhangAngle);
        // according to the paper the printing bed normal should be used here, but for xyz the printing bed normal is just the (in this case negative) printing direction (z-Axis)
        vtkVector3d printDir = lVtkHelper::invertVector(meshwa->printDir);
        
        // minimum z height needed to not be considered on bed
        double zTol = .000001;

        double ovhArea = 0;

        auto cellNorms = meshwa->getMeshCellNormals();

        // #pragma omp parallel for -> seems unsafe 
        for (vtkIdType i = 0; i < nCells; i++) {
            vtkVector3d cellNorm = vtkVector3d(cellNorms->GetTuple(i));
            auto curCosAngl = lVtkHelper::cosAngleNormVectors(cellNorm, printDir);

            // TODO: this is from paper, but really for 0? means ground plate?
            if (0 <= curCosAngl && curCosAngl <= criticalCosAngle) {
                // check if ground level / is on bed?
                vtkNew<vtkGenericCell> cll;
                poly->GetCell(i, cll);
                auto pts = cll->GetPoints();
                bool onBed = true;
                // is any point of the polygon above the bed level? If yes -> add area to cost
                for (vtkIdType ptIdx = 0; ptIdx < pts->GetNumberOfPoints(); ptIdx ++) {
                    double pt[3];
                    pts->GetPoint(ptIdx, pt);
                    if (abs(pt[2]) > zTol) {
                        onBed = false;
                        break;
                    }
                }
                if (onBed == false) {
                    // this ones 
                    auto triAr = cellAreas->GetComponent(i, 0);
                    ovhArea += triAr * fabs(curCosAngl);
                }
            } // if the facets normal is not in a critical angle (to flat overhang), nothing to do
        }

        if (printResults) std::cout << __func__ << " run with points: " << ovhArea << std::endl;
        return ovhArea;
    }

    // can easily be overflown as for a a flat angle / lots of considered as sharp angles
    // or a geometry with many tris the len of this tris adds up to a lot (a lot) of total length
    // (like intestines :D)
    double sharpness(MeshWithAttributes* meshwa) {

        // in a valid mesh at most two facets share a edge -> calulatinng something with respect to every neighbour of a facet
        // results in 2* the result as if on would only calculate it per edge -> /2 = per edge
        //auto cIt = polydata->NewCellIterator();
        //for (cIt->InitTraversal(); !cIt->IsDoneWithTraversal(); cIt->GoToNextCell()) {
        // upper bound in rad = 90° -> cehck paper for used value this ones the mathematical definition of sharp
        double sharpCornersAngle = meshwa->sharpCornersAngle;
        auto poly = meshwa->mesh;
        auto cellNorms = meshwa->getMeshCellNormals();
        vtkIdType nCells = poly->GetNumberOfCells();
        double cost = 0;
        // needed for the celledgeneighbors
        meshwa->createLinksIfNeeded();
        // -1: less than 3 edges, else last occurence of not == 1 neigbor or 1 if 1 neighbor
        // used for debug and visual, might be removed
        if (poly->GetCellData()->HasArray("NotOneNeighbor")) poly->GetCellData()->RemoveArray("NotOneNeighbor");
        vtkNew<vtkIdTypeArray> arNotOneNeighbor;
        arNotOneNeighbor->SetName("NotOneNeighbor");
        arNotOneNeighbor->SetNumberOfComponents(1);
        arNotOneNeighbor->SetNumberOfTuples(poly->GetNumberOfCells());
        // TODO: most likely safe
        //#pragma omp parallel for
        for (vtkIdType i = 0; i < nCells; i++) {
            
            auto curNorm = vtkVector3d(cellNorms->GetTuple(i));
            auto curCell = poly->GetCell(i);
            auto nEdges = curCell->GetNumberOfEdges();
            if (nEdges < 3) {
                arNotOneNeighbor->SetComponent(i, 0, -2);
                continue;
            }
            vtkIdType not3NNeighbors = 1;
            //assert(nEdges == 3); -> no because clipClosed may produce invalid, also may contain lines, points, ...
            for (vtkIdType k = 0; k < nEdges; k++) {
                auto curEdge = curCell->GetEdge(k);
                auto p0 = curEdge->GetPointIds()->GetId(0);
                auto p1 = curEdge->GetPointIds()->GetId(1);
                
                vtkNew<vtkIdList> cellNeighbors;
                poly->GetCellEdgeNeighbors(i, p0, p1, cellNeighbors);

                auto nNeigbors = cellNeighbors->GetNumberOfIds();
                if (nNeigbors != 1) {
                    not3NNeighbors = nNeigbors;
                }

                //assert(nNeigbors==1);
                auto nSharpNeig = 0;
                auto tempCost = 0;

                for (vtkIdType j = 0; j < nNeigbors; j++) {
                    auto neighborId = cellNeighbors->GetId(j);
                    auto neigNorm = vtkVector3d(cellNorms->GetTuple(neighborId));
                    auto curAngle = lVtkHelper::angleNormVectors(curNorm, neigNorm);
                    if (curAngle >= sharpCornersAngle) {
                        // sharp corner -> calculate the length of the edge and add it to the cost
                        double pc0 [3];
                        double pc1 [3];
                        // not threadsafe it says
                        poly->GetPoint(p0, pc0);
                        poly->GetPoint(p1, pc1);
                        double subVal[3];
                        vtkMath::Subtract(pc0, pc1, subVal);
                        auto edgeLen = vtkMath::Norm(subVal);
                        tempCost += edgeLen;
                        nSharpNeig += 1;
                    }
                }
                // remove effects that may occur if more than one neighbor exist
                // sadly this is what clipClosedSurfaceProduces
                if (nSharpNeig > 0) cost += tempCost / nSharpNeig;
            }
            //arNotOneNeighbor->InsertComponent(i, 0, not3NNeighbors);
            arNotOneNeighbor->SetComponent(i, 0, not3NNeighbors);
        }
        // used for debug and visual, might be removed
        poly->GetCellData()->AddArray(arNotOneNeighbor);
        // in a well defined tri mesh at most 2 facets share a edge.
        // As the cost is calulated for all facets for their neighbors every edge is considered twice
        // Therefore a edge that is sharp is added twice to the cost -> effect. doubling calc. cost.
        cost /= 2;

        if (printResults) std::cout << __func__ << " run with points: " << cost << std::endl;
        return cost;
    }

    double gap(MeshWithAttributes* meshwa) {
        // can all be replaced by a getter that calculates it if necessary else returns it
        // prequisites
        auto cellAreas = meshwa->getMeshTriAreas();
        auto gapCells = meshwa->getMeshGapCells();
        auto nGapCells = gapCells->GetNumberOfIds();

        vtkNew<vtkFloatArray> arGapCls;
        arGapCls->SetName("GapCellsCost");
        arGapCls->SetNumberOfComponents(1);
        arGapCls->SetNumberOfTuples(meshwa->mesh->GetNumberOfCells());
        arGapCls->Fill(-1);
        
        double cost = 0;
        //#pragma omp parallel for -> contains not threadsafe code as it seems (and might not even be worth it)
        for (vtkIdType i = 0; i < nGapCells; i++) {
            cost +=  cellAreas->GetComponent(gapCells->GetId(i), 0);
            arGapCls->SetComponent(gapCells->GetId(i), 0, cellAreas->GetComponent(gapCells->GetId(i), 0));
        }
        meshwa->mesh->GetCellData()->AddArray(arGapCls);

        if (printResults) std::cout << __func__ << " run with points: " << cost << std::endl;
        return cost;
    }

    double concavity(MeshWithAttributes* meshwa) {
        vtkNew<vtkMultiObjectMassProperties> maProp;
        auto poly = meshwa->mesh;
        auto ch = meshwa->convexHull;
        maProp->SetInputData(poly);
        maProp->Update();
        auto volPoly = maProp->GetTotalVolume();
        vtkNew<vtkMassProperties> maProp1;
        maProp1->SetInputData(ch);
        maProp1->Update();
        auto volCh = maProp1->GetVolume();
        auto conc = volCh - volPoly;
        if (printResults) std::cout << __func__ << " run with points: " << conc << std::endl;
        return conc;
    }

    // seems ok
    double feasibility(MeshWithAttributes* meshwa) {
        return meshwa->checkFeasibility();
    }

    // seems ok
    double feasibility(vtkSmartPointer<vtkPolyData> poly, double maxPrintDiagonal, double costFeasability) {
        auto diag = bounds2DiagLen(poly->GetBounds());
        if (diag > maxPrintDiagonal) return costFeasability;
        return 0;
    }

    double interfaceArea(MeshWithAttributes* meshwa) {
        auto poly = meshwa->mesh;
        // interface labels == 1 -> cell is on interface area / is new after cut
        auto ifLabels = poly->GetCellData()->GetArray("Labels");
        if (!ifLabels) {
            //n TODO: check what to do
            cout << "interfaceArea needs data with Labels array. (vtkClipClosedSurface->SetScalarModeToLabels())" << endl;
            return 0;
            // throw "interfaceArea needs data with Labels array. (vtkClipClosedSurface->SetScalarModeToLabels())";
        }
        auto nTups = ifLabels->GetNumberOfTuples();
        double ifArea = 0;

        for (vtkIdType i = 0; i < nTups; i++) {
            // is one for interface cells
            if (ifLabels->GetComponent(i, 0) > .5) {
                vtkNew<vtkGenericCell> tempCell;
                poly->GetCell(i, tempCell);
                ifArea += lVtkHelper::computeTriArea(tempCell);
            }
        }
        if (printResults) std::cout << __func__ << " run with points: " << ifArea << std::endl;
        return ifArea;
    }

    // calculates the number of unconnected regions in the part / 
    // for good meshes the number of parts
    double quantity(MeshWithAttributes* meshwa) {
        // TODO:
        // there were cases where the mesh looks good after no cut but is fully non manifold resultinh in every facet is a part
        // if it is fixable by clean polydata it is done at the output and therefore ok to do here for better results
        // one should probably find out where this non manifoldness results from -> todo
        vtkNew<vtkPolyDataConnectivityFilter> connFilter;
        connFilter->SetInputData(meshwa->mesh);
        connFilter->SetExtractionModeToAllRegions();
        connFilter->Update();
        double nRegions = static_cast<double>(connFilter->GetNumberOfExtractedRegions());
        if (printResults) std::cout << __func__ << " run with points: " << nRegions << std::endl;
        return nRegions;
    }


    //---------------------------------------------------------------------------------------------------
    /* #endregion quality measures as desctibed ... */

}