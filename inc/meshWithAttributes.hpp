/**
 * @file meshWithAttributes.hpp
 * @author your name (you@domain.com)
 * @brief a header only class as data container with some additional functionality. May contain mesh and convex hull(CH)
 * @version 0.1
 * @date 2022-03-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef MESH_WITH_ATTRIBUTES_HPP
#define MESH_WITH_ATTRIBUTES_HPP

#include <vector>
#include <set>
// needed for windows to use constants
#define _USE_MATH_DEFINES
#include <math.h>
#include <thread>
// for timing tests
#include <chrono>
#include <limits>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkVector.h>
#include <vtkOBBTree.h>
#include <vtkDataArray.h>
#include <vtkCellData.h>
#include <vtkGenericCell.h>
#include <vtkPlane.h>
#include <vtkPlaneCollection.h>
#include <vtkModifiedBSPTree.h>
#include <vtkAbstractCellLocator.h>
#include <vtkStaticCellLocator.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkCellTreeLocator.h>
#include <vtkPoints.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkClipClosedSurface.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>
#include <vtkCenterOfMass.h>
#include <vtkBoundingBox.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>

#include "typesConstantsDefinition.hpp"
#include "lVtkHelper.hpp"
#include "arrayMath.hpp"

#include "SimpleLogger.hpp"

namespace polyQualityMeasures
{

    class MeshWithAttributes
    {
        // most public, private is only needed in special cases and makes modifiing / composition harder
    public:
        // the mesh / part
        vtkSmartPointer<vtkPolyData> mesh;

        // the printing direction of the mesh (z-Axis for normal 3 Axis xyz)
        // defaults to upright print dir (0, 0, 1)
        vtkVector3d printDir = vtkVector3d(0, 0, 1);

        // overhang threshold angle in rad (default 45* M_PI / 180)
        double overhangAngle = 45 * M_PI / 180;

        // for convex hulls to find the printing direction and the greatest facet
        vtkIdType greatestCellId;
        double greatestCellSize = 0;

        // the cost for unprintable bodies (diag to big)
        double costFeasability = 1000;

        // this is a per printer attribute defining if  the printer could theoretically print the obj in mm³
        // -> ! this is only a upper printability limit -> above the part can definitly not be printed.
        //      below it might be possible to print it, might also not be possible (one len > axis ranges)
        double maxPrintDiagonal = 100 * 100 * 100;

        // the printers dimension -> currently only works with box like build volumes.
        // there exist also cylinder like and other types of build volumes. these would need to be implemented if needed
        vtkBoundingBox printerBuildVolume;

        // the angle below which angles are considered sharp in rad (from paper)
        double sharpCornersAngle = M_PI * (170. / 180.);

        // the maxiumum length of the normal from the center of a face before hitting an other face to be considered a non-allowable gap
        // got to check if the paper proposed a value, but .5 is in the range of one nozzle width -> should be easily printable
        // from paper 2 mm which semms a bit high but ...
        double gapDistance = 2;

        // some operations need vtk buildLinks executed on the mesh
        // this variable contains if the links have already been build
        bool linksCreated = false;

        // some tolerance, where used, values below this are considered as equal
        double lowerTol = 1e-13;

        vtkSmartPointer<vtkIdList> gapCells;
        std::vector<std::vector<double>> vecGapCellBnds;

        // the mesh / part
        // the convex hull of the part if one is set
        vtkSmartPointer<vtkPolyData> convexHull;
        // the cell id of the greatest cell of the convex hull if calcualted
        vtkIdType greatestCHCellId;
        // the size of the greatest cell of the ch if calculated
        double greatestCHCellSize = 0;
        // multiple trees are defined because one can change the used one, could prbly be replaced by some parent class ptr
        // the obb tree for gap / collision detection / ray casting
        vtkSmartPointer<vtkOBBTree> obbTree;
        // tree for gap / collision detection / ray casting
        vtkSmartPointer<vtkModifiedBSPTree> modifiedBspTree;
        // tree for gap / collision detection / ray casting
        vtkSmartPointer<vtkCellTreeLocator> cellTreeLocator;
        // tree for gap / collision detection / ray casting
        vtkSmartPointer<vtkStaticCellLocator> staticCellLocator;

        vtkSmartPointer<vtkCellTreeLocator> chCellTreeLocator;

        MeshWithAttributes()
        {
            vtkBoundingBox bv;
            bv.SetBounds(0, 100, 0, 100, 0, 100);
            this->setBuildVolume(bv);
        }

        void setBuildVolume(vtkBoundingBox bv)
        {
            this->maxPrintDiagonal = bv.GetDiagonalLength();
            this->printerBuildVolume = bv;
        }

        void setBuildVolume(double xLen, double yLen, double zLen)
        {
            vtkBoundingBox bv;
            bv.SetBounds(0, xLen, 0, yLen, 0, zLen);
            this->setBuildVolume(bv);
        }

        struct sGapPtsAndCost
        {
            // contains startP ans intersectP alternating enabling to check if n==even and n+1==n+1 -> intersection in bounds
            vtkSmartPointer<vtkPoints> gapStartAndIntersectPts;
            std::vector<double> gapCosts;
        };

        double gapCost = -1;

        // information about the geps that were found after compute gap tris
        struct sGapPtsAndCost gapPtsAndCost;

        /* #region getter's */

        // I (RULE OF THREE)
        ~MeshWithAttributes()
        {
            ;
        }

        // copy constructor II
        MeshWithAttributes(const MeshWithAttributes &mwa)
        {
            this->copyMwa(mwa);
        }

        // III
        MeshWithAttributes &operator=(const MeshWithAttributes &original)
        {
            if (this == &original)
                return *this;

            copyMwa(original, *this);

            return *this;
        }

        MeshWithAttributes &copyMwa(const MeshWithAttributes &original, MeshWithAttributes &target)
        {
            // "allocate"
            auto newMesh = vtkSmartPointer<vtkPolyData>::New();
            auto newCh = vtkSmartPointer<vtkPolyData>::New();
            // "populate"
            // reallocate
            if (original.mesh)
            {
                newMesh->DeepCopy(original.mesh);
                target.mesh = newMesh;
            }
            if (original.convexHull)
            {
                newCh->DeepCopy(original.convexHull);
                target.convexHull = newCh;
            }

            // some could be bundeled :| -> here copyWithoutMeshes should be called but no risk for now :D
            target.printDir = original.printDir;
            target.overhangAngle = original.overhangAngle;
            target.greatestCellId = original.greatestCellId;
            target.greatestCellSize = original.greatestCellSize;
            target.costFeasability = original.costFeasability;
            target.maxPrintDiagonal = original.maxPrintDiagonal;
            target.sharpCornersAngle = original.sharpCornersAngle;
            target.gapDistance = original.gapDistance;
            target.linksCreated = original.linksCreated;
            target.greatestCHCellId = original.greatestCHCellId;
            target.greatestCHCellSize = original.greatestCHCellSize;
            target.setBuildVolume(original.printerBuildVolume);
            /*std::cout << "orig clls: " << original.mesh->GetNumberOfCells() << std::endl;
            std::cout << "targ clls: " << target.mesh->GetNumberOfCells() << std::endl;
            std::cout << "cpy dn" << std::endl;*/
            return target;
        }

        MeshWithAttributes &copyMwa(const MeshWithAttributes &original)
        {
            return copyMwa(original, *this);
        }

        // TODO: not sure if all is in -> this is erroneous causing segfs when used
        /*    MeshWithAttributes operator=(const MeshWithAttributes& original) {
                MeshWithAttributes cloneObj(original);
                //cloneObj.copyMwa(original);
                return cloneObj;
            }*/

        MeshWithAttributes &copyWithoutMeshes(const MeshWithAttributes &original, MeshWithAttributes &target)
        {
            // some could be bundeled :|
            target.printDir = original.printDir;
            target.overhangAngle = original.overhangAngle;
            target.greatestCellId = original.greatestCellId;
            target.greatestCellSize = original.greatestCellSize;
            target.costFeasability = original.costFeasability;
            target.maxPrintDiagonal = original.maxPrintDiagonal;
            target.sharpCornersAngle = original.sharpCornersAngle;
            target.gapDistance = original.gapDistance;
            target.setBuildVolume(original.printerBuildVolume);

            return target;
        }

        MeshWithAttributes &copyWithoutMeshes(const MeshWithAttributes &original)
        {
            return copyWithoutMeshes(original, *this);
        }

        /**
         * @brief Get the To Print Dir Rotated Mesh object
         *
         * @param appliedTransform should assing the input pointer to the used transform if this is correct and it does not go out of scope
         * @return vtkSmartPointer<vtkPolyData>
         */
        vtkSmartPointer<vtkPolyData> getToPrintDirRotatedMesh(vtkSmartPointer<vtkTransform> appliedTransform = nullptr)
        {
            auto prntTrans = this->getPrintDirAsTransform();
            vtkNew<vtkTransformPolyDataFilter> rotFltr;
            rotFltr->SetInputData(this->mesh);
            rotFltr->SetTransform(prntTrans);
            rotFltr->Update();
            vtkSmartPointer<vtkPolyData> rotPoly = rotFltr->GetOutput();
            if (appliedTransform != nullptr)
            {
                appliedTransform->SetMatrix(prntTrans->GetMatrix());
            }
            return rotPoly;
        }

        /**
         * @brief Get the Bbx Aligned Mesh object. Align the poly bbx with the printer build volume by rotation around z axis, so that the long x or y side of the poly is aligned with the long x or y side of the printer
         *
         * @param poly
         * @return vtkSmartPointer<vtkPolyData>
         */
        vtkSmartPointer<vtkPolyData> getBbxAlignedMesh(vtkSmartPointer<vtkPolyData> poly, vtkSmartPointer<vtkTransform> appliedTransform = nullptr)
        {
            auto zTransf = getAlignBbxsTransform(this->printerBuildVolume, getPolyBoundsAsVtkBoundingBox(poly));
            vtkNew<vtkTransformPolyDataFilter> rotZFltr;
            rotZFltr->SetInputData(poly);
            rotZFltr->SetTransform(zTransf);
            rotZFltr->Update();
            auto alignPoly = vtkSmartPointer<vtkPolyData>::New();
            alignPoly = rotZFltr->GetOutput();
            if (appliedTransform != nullptr)
            {
                appliedTransform->SetMatrix(zTransf->GetMatrix());
            }
            // TODO: remove debug only
            return alignPoly;
        }

        /**
         * @brief Get the transform (0, 90°) that alignes the longest x or y side of the printer with the parts longest x or y side
         *
         * @param printerBbx the bbx of the printer
         * @param partBbx the bbx of the part
         * @return vtkSmartPointer<vtkTransform> a rotation transform
         */
        vtkSmartPointer<vtkTransform> getAlignBbxsTransform(vtkBoundingBox printerBbx, vtkBoundingBox partBbx)
        {
            auto aligAngle = lVtkHelper::getAlignBbxsAngle(printerBbx, partBbx);
            auto aligTransf = vtkSmartPointer<vtkTransform>::New();
            aligTransf->RotateZ(aligAngle);
            return aligTransf;
        }

        /**
         * @brief Get the Poly Bounds As Vtk Bounding Box object
         *
         * @param poly the mesh
         * @return vtkBoundingBox the bbx of poly
         */
        vtkBoundingBox getPolyBoundsAsVtkBoundingBox(vtkSmartPointer<vtkPolyData> poly)
        {
            double bnds[6];
            poly->GetBounds(bnds);
            vtkBoundingBox bbx;
            bbx.SetBounds(bnds);
            return bbx;
        }

        /**
         * @brief Get the Origin Translation For Rotated Mesh object, meaning:
         *  after the mesh got rotated so the print direction is the z direction ("rotate so that print direction vector and z direction vector are the same")
         *  the mesh needs to be translated so it's start coordinates are (0, 0, 0) and end coordinates are (xWidth, yWitdh, zWidth)
         *  this is needed so the part can be printed easily without any further steps also i think there is something in the 3mf specs.
         *
         * @param bool needsRotation: if the part needs to be rotated to z direction before calculation the translation or not. (default=true, because less harmfull if wrong)
         *
         * @return vtkTransform: the transform of the translation operation needed to bring the part to the origin (0,0,0)
         */
        vtkSmartPointer<vtkTransform> getOriginTranslationForRotatedMesh(vtkSmartPointer<vtkPolyData> poly, bool needsRotation = true, bool alignBbxs = true, vtkSmartPointer<vtkTransform> appliedRotationTransform = nullptr, vtkSmartPointer<vtkTransform> appliedAlignBbxsTransform = nullptr, vtkSmartPointer<vtkTransform> totalTransform = nullptr)
        {
            // TODO
            // create the pointers always so they exist if needed
            if (appliedRotationTransform == nullptr)
                appliedRotationTransform = vtkSmartPointer<vtkTransform>::New();
            if (appliedAlignBbxsTransform == nullptr)
                appliedAlignBbxsTransform = vtkSmartPointer<vtkTransform>::New();

            auto oriPoly = vtkSmartPointer<vtkPolyData>::New();
            if (needsRotation)
            {
                oriPoly = getToPrintDirRotatedMesh(appliedRotationTransform);
            }
            else
            {
                oriPoly = poly;
            }

            auto aliPoly = vtkSmartPointer<vtkPolyData>::New();
            if (alignBbxs)
            {
                aliPoly = getBbxAlignedMesh(oriPoly, appliedAlignBbxsTransform);
            }
            else
            {
                aliPoly = oriPoly;
            }

            double bnds[6];
            aliPoly->GetBounds(bnds);
            vtkNew<vtkTransform> translTransf;
            // lower toll is needed to (mostly) ensure that the lower part is at 0 or on the positive side of 0 near 0
            translTransf->Translate(-bnds[0] + lowerTol, -bnds[2] + lowerTol, -bnds[4] + lowerTol);
            vtkNew<vtkTransformPolyDataFilter> transfFltr;
            // std::cout << "transl:" << std::endl;
            // lVtkHelper::printStringFromvtk4x4MatrixValues(translTransf);

            if (totalTransform != nullptr)
            {
                // assuming premultiply
                totalTransform->Concatenate(translTransf);
                totalTransform->Concatenate(appliedAlignBbxsTransform);
                totalTransform->Concatenate(appliedRotationTransform);
            }

            return translTransf;
        }

        vtkSmartPointer<vtkTransform> getOriginTranslationForRotatedMesh()
        {
            return getOriginTranslationForRotatedMesh(this->mesh, true);
        }

        vtkSmartPointer<vtkTransform> getOriginTranslationForRotatedMesh(vtkTransform *totalTransform)
        {
            return getOriginTranslationForRotatedMesh(this->mesh, true, true, nullptr, nullptr, totalTransform);
        }

        vtkSmartPointer<vtkTransform> getPrintRotationAndOriginTranslationAsTransform()
        {
            auto totalTransf = vtkSmartPointer<vtkTransform>::New();
            // auto printDirRot = this->getPrintDirAsTransform();
            // auto printDirTransl = getOriginTranslationForRotatedMesh(totalTransf.GetPointer());
            // printDirRot->Concatenate(printDirTransl);
            // return printDirRot;
            return totalTransf;
        }

        /**
         * @brief Get the Mesh Origined Rotated To Print Dir object     orients this mesh so that the print dir is the global z dir and the longest side of the printer and the part are aligned
         *
         * totalApplied... -> getter for applied transf
         *
         * @return vtkSmartPointer<vtkPolyData>
         */
        vtkSmartPointer<vtkPolyData> getMeshOriginedRotatedToPrintDir(vtkSmartPointer<vtkTransform> totalAppliedTransform = nullptr)
        {
            auto rotTransf = vtkSmartPointer<vtkTransform>::New();
            auto rotPoly = getToPrintDirRotatedMesh(rotTransf);

            // total could also be used instead of this get 2 and conc
            auto bbxAlignTransf = vtkSmartPointer<vtkTransform>::New();
            auto translTransf = getOriginTranslationForRotatedMesh(rotPoly, false, true, nullptr, bbxAlignTransf);
            translTransf->Concatenate(bbxAlignTransf);
            vtkNew<vtkTransformPolyDataFilter> transfFltr;
            transfFltr->SetInputData(rotPoly);
            transfFltr->SetTransform(translTransf);
            transfFltr->Update();
            vtkSmartPointer<vtkPolyData> orientPoly = transfFltr->GetOutput();

            if (totalAppliedTransform != nullptr)
            {
                totalAppliedTransform->Concatenate(translTransf);
                totalAppliedTransform->Concatenate(rotTransf);
            }

            return orientPoly;
        }

        /**
         * @brief get a cell locator. If it does not exist, create it, also creates links if needed
         *
         * @tparam T the locator type
         * @param targetPointer the pointer to which to assign the created object
         * @param force if the creation / recreation is requested even if it already exists
         * @return vtkSmartPointer<T> the pointer to the locator
         */
        template <typename T>
        vtkSmartPointer<T> safeGetLocator(vtkSmartPointer<T> targetPointer, bool force = false)
        {
            if (targetPointer == nullptr || force)
            {
                this->createLinksIfNeeded();
                auto ot = vtkSmartPointer<T>::New();
                ot->SetDataSet(this->mesh);
                targetPointer = ot;
                targetPointer->BuildLocator();
            }
            return targetPointer;
        }

        // safe access, tree gets created if neccessary
        vtkSmartPointer<vtkOBBTree> getObbTree(bool force = false)
        {
            return safeGetLocator(this->obbTree, force);
        }

        // safe access, tree gets created if neccessary
        vtkSmartPointer<vtkModifiedBSPTree> getModifiedBspTree(bool force = false)
        {
            return safeGetLocator(this->modifiedBspTree, force);
        }

        // safe access, tree gets created if neccessary
        vtkSmartPointer<vtkCellTreeLocator> getCellTreeLocator(bool force = false)
        {
            return safeGetLocator(this->cellTreeLocator, force);
        }

        // safe access, tree gets created if neccessary
        vtkSmartPointer<vtkStaticCellLocator> getStaticCellLocator(bool force = false)
        {
            return safeGetLocator(this->staticCellLocator, force);
        }

        // enables easy change of default locator
        // for the test model vtkModifiedBspTree is a lot faster than vtkObbTree (x3),
        // changed to static cell locator, because of segfs with the other locator
        // seems to be as fast or faster
        // works only if the same interfaces are supported
        vtkSmartPointer<vtkStaticCellLocator> getCellLocator(bool force = false)
        {
            return this->getStaticCellLocator(force);
        }

        /**
         * @brief safe access to the mesh normals of this mesh, if none exist they get calulated
         *
         * @param force if normals should be created even if they already exist
         * @return vtkSmartPointer<vtkDataArray> the normals as array
         */
        vtkSmartPointer<vtkDataArray> getMeshCellNormals(bool force = false)
        {
            if (this->mesh->GetCellData()->GetNormals() == nullptr || force)
            {
                lVtkHelper::calcualteCellNormalsInPlace(this->mesh);
            }
            return this->mesh->GetCellData()->GetNormals();
        }

        // safe access to the ch normals, if none exist they get calulated
        vtkSmartPointer<vtkDataArray> getChCellNormals()
        {
            if (this->convexHull->GetCellData()->GetNormals() == nullptr)
            {
                lVtkHelper::calcualteCellNormalsInPlace(this->convexHull);
            }
            return this->convexHull->GetCellData()->GetNormals();
        }

        // calculates all tri areas if no array with name defined in TRI_AREA_ARR_NAME exists (typesConstantsDefinition)
        vtkDataArray *getMeshTriAreas(bool force = false)
        {
            auto hasCellAreas = this->mesh->GetCellData()->HasArray(TRI_AREA_ARR_NAME);
            if (hasCellAreas == 0 || force)
            {
                lVtkHelper::appendCellAreasArrayIp(this->mesh);
            }
            return this->mesh->GetCellData()->GetArray(TRI_AREA_ARR_NAME);
        }

        /**
         * @brief calculates the cell centers if needed, they get added as an cell array with name defined in CELL_CENTER_ARR_NAME
         *
         * @param force if creation is forced even if they already exist
         * @return vtkDataArray* the cell center data
         */
        vtkDataArray *getMeshCellCenters(bool force = false)
        {
            auto cellCenters = this->mesh->GetCellData()->HasArray(CELL_CENTER_ARR_NAME);
            if (cellCenters == 0 || force)
            {
                this->mesh = lVtkHelper::appendCellCenters(this->mesh);
            }
            return this->mesh->GetCellData()->GetArray(CELL_CENTER_ARR_NAME);
        }

        /**
         * @brief get the ids of the cells that have a gap which is to small according to the gap criteria (normal ray from center -> collision = distance)
         *
         * @param force
         * @return vtkSmartPointer<vtkIdList>
         */
        vtkSmartPointer<vtkIdList> getMeshGapCells(bool force = false)
        {
            if (this->gapCells == nullptr || force)
            {
                this->computeGapTris();
            }
            return this->gapCells;
        }

        /**
         * @brief get the id of the gretest cell of the ch. Calculate if needed
         *
         * @return vtkIdType
         */
        vtkIdType getGreatestCHCellId()
        {
            if (greatestCHCellSize == 0)
                computeCHGreatestCellId();
            return this->greatestCHCellId;
        }

        /**
         * @brief check if 2 planes are equal up to some errorMargin, equal: same normal and points are one one plane with that normal
         *
         * @param x first plane
         * @param y second plane
         * @return true are equal
         * @return false are not equal
         */
        bool planesEqual(const vtkSmartPointer<vtkPlane> &x,
                         const vtkSmartPointer<vtkPlane> &y) const
        {
            // "()" avoid macro expansion on windows which causes failed compilation
            auto errorMargin = (std::numeric_limits<double>::min)() * 10;
            // auto errorMargin = 1;
            //  point on plane
            auto px = x->GetOrigin();
            auto py = y->GetOrigin();
            // plane normal
            auto nx = x->GetNormal();
            auto ny = y->GetNormal();
            //  size is 3 for pts and normals
            for (auto i = 0; i < 3; i++)
            {
                auto tnx = nx[i];
                auto tny = ny[i];
                // check if similar
                if (fabs(tnx - tny) > errorMargin)
                {
                    return false;
                }
            }
            // if the code ran until here, the normal is similar -> check if the points lie most likely on the same plane
            // in same plane: dot(A-B, n) = 0 (a or b needs to lie on plane) -> cosangle is 0 -> angle is n*90deg
            // check if point y lies in x
            for (auto i = 0; i < 3; i++)
            {
                auto dp = vtkPlane::Evaluate(nx, px, py);
                if (fabs(dp) > errorMargin)
                    return false;
            }
            return true;
        }

        /**
         * @brief check if planeVec contains a plane that is equal to new plane, equal: same normal and point on the same plane
         *
         * @param planeVec
         * @param newPlane
         * @return true
         * @return false
         */
        bool planeInVec(const std::vector<vtkSmartPointer<vtkPlane>> &planeVec, const vtkSmartPointer<vtkPlane> newPlane)
        {
            for (auto pln : planeVec)
            {
                if (planesEqual(pln, newPlane))
                    return true;
            }
            return false;
        }

        /**
         * @brief insert new plane to plane vector only if no equal plane exists, equal: same normal and point on the same plane
         *
         * @param planeVec
         * @param newPlane
         * @return true
         * @return false
         */
        bool insertUniquePlaneInVec(std::vector<vtkSmartPointer<vtkPlane>> &planeVec, vtkSmartPointer<vtkPlane> newPlane)
        {
            if (planeInVec(planeVec, newPlane))
                return false;
            planeVec.push_back(newPlane);
            return true;
        }

        /**
         * @brief get the convex hull as set of unique planes / bounding planes
         *
         * @param invertNormals the normal direction. If they should point inwards (prbly true) or outwards
         * @param filterOutNoCutPlanes if the resulting planes should be filtered so only unique planes remain
         * @return vtkSmartPointer<vtkPlaneCollection>
         */
        vtkSmartPointer<vtkPlaneCollection> getChCellsAsPlanes(bool invertNormals = true, bool filterOutNoCutPlanes = true)
        {
            auto nCells = this->convexHull->GetNumberOfCells();
            auto norms = this->getChCellNormals();

            vtkNew<vtkPlaneCollection> planeCol;
            vtkNew<vtkPlaneCollection> planeCol1;
            double p[3];
            double n[3];

            // refac to either struct or just vect of tuple or vect[6]
            std::vector<vtkSmartPointer<vtkPlane>> planeSet;

            vtkNew<vtkGenericCell> cell;
            for (vtkIdType i = 0; i < nCells; i++)
            {
                this->convexHull->GetCell(i, cell);
                cell->GetPoints()->GetPoint(0, p);
                norms->GetTuple(i, n);
                vtkNew<vtkPlane> plane;
                plane->SetOrigin(p);
                if (invertNormals)
                    mul3Ip(n, -1.);
                plane->SetNormal(n);
                planeCol->AddItem(plane);
                // struct sPlane pln;
                // pln.origin = plane->GetOrigin();
                // pln.normal = plane->GetNormal();
                insertUniquePlaneInVec(planeSet, plane);
                // auto planeTup = std::make_tuple(p[0], p[1], p[2], n[0], n[1], n[2]);
                // planeSet.insert(planeTup);
            }

            // create unique set of planes
            for (auto pln : planeSet)
            {
                planeCol1->AddItem(pln);
            }

            if (filterOutNoCutPlanes)
            {
                vtkNew<vtkPlaneCollection> fltrPlanc;
                lVtkHelper::filterPlaneCol4Cut(this->mesh->GetPoints()->GetData(), planeCol1, fltrPlanc);
                return fltrPlanc;
            }

            return planeCol1;
        }

        /* #endregion getter's */

        /**
         * @brief Create a Links If Needed for the mesh, f.e. needed for cell locator
         *
         * @param force
         * @return true
         * @return false
         */
        bool createLinksIfNeeded(bool force = false)
        {
            if (!this->linksCreated || force)
            {
                this->mesh->BuildLinks();
                this->linksCreated = true;
                return true;
            }
            return false;
        }

        /**
         * @brief cut this->mesh by this->convexHull ! might produce not watertight bodys
         *
         */
        void cutMeshWithCh()
        {
            // TODO
            auto clipPlanes = this->getChCellsAsPlanes();

            vtkNew<vtkPlaneCollection> newPlaneCol;
            lVtkHelper::filterPlaneCol4Cut(this->mesh->GetPoints()->GetData(), clipPlanes, newPlaneCol);
            if (newPlaneCol->GetNumberOfItems() == 0)
            {
                return;
            }

            vtkNew<vtkClipClosedSurface> clipClosedSurf;
            clipClosedSurf->SetInputData(this->mesh);
            // TODO TriangulationErrorDisplay
            clipClosedSurf->SetClippingPlanes(newPlaneCol);
            // e options are "None", "Colors", and "Labels". For the "Labels" option, a scalar value of "0" indicates an original cell, "1" indicates a new cell on a cut face, and "2" indicates a new cell on the ActivePlane as set by the SetActivePlane() method. The default scalar mode is "None".
            clipClosedSurf->SetScalarModeToLabels();
            clipClosedSurf->PassPointDataOn();
            clipClosedSurf->Update();
            clipClosedSurf->SetTriangulationErrorDisplay(true);
            auto triangError = clipClosedSurf->GetErrorCode();
            if (triangError)
            {
                logger.addErrLog("clip closed failed in creating closed");
                logger.printLastErr();
            }
            auto cleanedMesh = lVtkHelper::getCleanPolydata(clipClosedSurf->GetOutput());
            this->mesh = cleanedMesh;

            /*
            //if this is used mallocs get thrown
            auto clpMsh = vtkSmartPointer<vtkPolyData>::New();
            clpMsh = lVtkHelper::clipClosedKeepCellLabels(this->mesh, clipPlanes);
            //auto cleanMesh = vtkSmartPointer<vtkPolyData>::New();
            //cleanMesh = lVtkHelper::getCleanPolydata(clpMsh);
            this->mesh = clpMsh;
            //this->mesh = cleanMesh;*/
        }

        /**
         * @brief Get the Print Dir As Transform object
         *
         * @return vtkSmartPointer<vtkTransform>
         */
        vtkSmartPointer<vtkTransform> getPrintDirAsTransform()
        {
            // the aim is to get a transform, so that the print dir is the new z-axis
            // this is done by calculating the angles in relation to the x and z-axis.
            // inverting both transforms the vector will be returned to the z-axis
            auto pd = this->printDir;
            pd.Normalize();
            // relative to the xAxis (around z) (math. positive in relation to the x axis)
            double xAngle;
            if (pd.GetY() > (1 - (std::numeric_limits<double>::min)() * 100) && fabs(pd.GetX()) < ((std::numeric_limits<double>::min)() * 100))
            {
                // only y no x case
                // degree use only positive anlges -> -90 = 270°
                xAngle = 90 * (pd.GetY() > 0 ? 1 : 3);
            }
            else if (pd.GetY() < (1 - (std::numeric_limits<double>::min)() * 100))
            {
                // no y case check if no x also
                xAngle = (fabs(pd.GetX()) < ((std::numeric_limits<double>::min)() * 100)) ? 0. : 90.;
            }
            else
            {
                // y is not 0
                xAngle = atan(pd.GetX() / pd.GetY());
                xAngle = xAngle / M_PI * 180;
            }
            // relative to the xy-plane (around the axis normal to the vec created after yzAngle) / new x
            // arround the new x axis starting from the z-axis (see wikipedia Kugelkoordinaten)
            auto zAngle = acos(pd.GetZ()) / M_PI * 180;
            vtkNew<vtkTransform> transform;
            transform->RotateZ(xAngle);
            // transform->PostMultiply();
            transform->RotateX(zAngle);
            return transform;
        }

        // cant use polyquality... -> include loop
        /**
         * @brief check if the poly is feasible / the diagonal is smaller than the max print diagonal
         *
         * @param poly
         * @param maxPrintDiagonal
         * @param costFeasability
         * @return double
         */
        double checkFeasibility(vtkSmartPointer<vtkPolyData> poly, double maxPrintDiagonal, double costFeasability)
        {
            auto diag = bounds2DiagLen(poly->GetBounds());
            if (diag > maxPrintDiagonal)
                return costFeasability;
            return 0;
        }

        // cant use polyquality... -> include loop
        double checkFeasibility(MeshWithAttributes *meshwa)
        {
            return this->checkFeasibility(this->mesh, this->maxPrintDiagonal, this->costFeasability);
        }

        double checkFeasibility(vtkSmartPointer<vtkPolyData> poly)
        {
            return this->checkFeasibility(poly, this->maxPrintDiagonal, this->costFeasability);
        }

        double checkFeasibility()
        {
            return this->checkFeasibility(this);
        }

        // not the most efficient way of calculating these things but easily readable. The compiler opts it prbly anyway
        /**
         * @brief calcualte the number of cuts, as real number, needed so that original len fits in target len
         *
         * @param originalLen
         * @param targetLen
         * @return double
         */
        double calculateRcuts(double originalLen, double targetLen)
        {
            // std::cout << "origLen: " << originalLen << " targetLen: " << targetLen << " prod: " << (originalLen/targetLen) << std::endl;
            return (originalLen / targetLen) - 1;
        };
        // wrapper for semantic reasons
        /**
         * @brief calcualte the integer number of cuts from the rcuts
         *
         * @param valRcuts
         * @return double
         */
        double calculateNcuts(double valRcuts) { return ceil(valRcuts); };
        /**
         * @brief calculate the integer number of cuts needed so that original len fits in target len
         *
         * @param originalLen
         * @param targetLen
         * @return double
         */
        double calculateNcuts(double originalLen, double targetLen) { return ceil(calculateRcuts(originalLen, targetLen)); };
        // remainder length after cutting the part (lenght of the last part if no offset is used)
        double calculateLRem(double originalLen, double targetLen)
        {
            double valRCuts = calculateRcuts(originalLen, targetLen);
            auto rem = valRCuts - floor(valRCuts);
            auto valRemLen = rem * targetLen;
            return valRemLen;
        };
        double calculateLOffset(double originalLen, double targetLen)
        {
            return calculateLRem(originalLen, targetLen) / 2;
        };

        /**
         * @brief check if inner bbx fits in outer bbx, with some upperTol erance added to the size of the outer bbx
         *
         * @param innerBbx
         * @param outerBbx
         * @return true
         * @return false
         */
        bool bbxFitsInBbx(vtkBoundingBox innerBbx, vtkBoundingBox outerBbx)
        {
            // check if inner bbx would fit in outer bbx uses upperTol
            return innerBbx.GetLength(0) < (outerBbx.GetLength(0) * upperTol) && innerBbx.GetLength(1) < (outerBbx.GetLength(1) * upperTol) && innerBbx.GetLength(2) < (outerBbx.GetLength(2) * upperTol);
        }

        double upperTol = .9999999;

        // cut along the direction unitl the poly fits in the target bbx with the giben direction (0:x, 1:y, 2:z)
        std::vector<vtkSmartPointer<vtkPolyData>> cutAlongAxisUntilItFits(vtkSmartPointer<vtkPolyData> poly, double originalLen, double targetLen, int dir)
        {
            // dir: x, y, z (0,1,2)
            // the part is to tall

            // i thoudht the cuts should always be at the same place, but that might create thin slices -> caluclate the original len dynamically for the part.
            // well no i did it so the ch will also fit..... the original len should already be part based
            targetLen *= this->upperTol;
            auto nCuts = calculateNcuts(originalLen, targetLen);
            // the order here is mixed up
            // auto nCuts = calculateNcuts(origLen, targetLen);

            // std::cout << "nCuts: " << nCuts << std::endl;
            double arDir[3] = {0, 0, 0};
            double negArDir[3] = {0, 0, 0};
            arDir[0] = (dir == 0) ? 1 : 0;
            arDir[1] = (dir == 1) ? 1 : 0;
            arDir[2] = (dir == 2) ? 1 : 0;
            negArDir[0] = (dir == 0) ? -1 : 0;
            negArDir[1] = (dir == 1) ? -1 : 0;
            negArDir[2] = (dir == 2) ? -1 : 0;

            // off also needs to be offset by the starting position of the mesh in the chosen direction
            double bnds[6];
            poly->GetBounds(bnds);
            // idx: 0, 2, 4
            double startOff = bnds[dir * 2];

            // start an offset before (rem/2) -> only if small last fill ( < .5)
            // auto rcuts = calculateRcuts(originalLen, targetLen);
            auto off = 0;
            /*if (rcuts < .4) {
                off = calculateLOffset(originalLen, targetLen);
            }*/
            off = calculateLOffset(originalLen, targetLen);
            if (fabs(off) < 1)
                off = 0;

            std::vector<vtkSmartPointer<vtkPolyData>> vecClipedPolys;

            // if nothing to do
            /*std::cout << "\n\nnCuts: " << nCuts << " orig: " << originalLen << " tarlen: " << targetLen << " off: " << off << " startOff: " << startOff << std::endl;
            std::cout << "bnds: " << array2Str6d(bnds) << " rcuts: " << rcuts << std::endl;
            std::cout << "\n" << std::endl;*/
            if (nCuts < 1.)
            {
                // std::cout << "return" << std::endl;
                vecClipedPolys.push_back(poly);
                return vecClipedPolys;
            }

            // auto off = calculateLOffset(origLen, targetLen);
            vtkSmartPointer<vtkPolyData> curPoly = poly;
            for (auto i = 0; i < nCuts; i++)
            {
                // do them cuts -> might create == cases (len = max len)
                auto cutPos = ((i + 1) * targetLen) - off + startOff;
                // std::cout << "cutPos: " << cutPos << " tarLen: " << targetLen << " off: " << off << " startOff: " << startOff << std::endl;
                double arCutPos[3] = {0, 0, 0};
                arCutPos[0] = (dir == 0) ? cutPos : 0;
                arCutPos[1] = (dir == 1) ? cutPos : 0;
                arCutPos[2] = (dir == 2) ? cutPos : 0;
                vtkNew<vtkPlaneCollection> plCol;
                vtkNew<vtkPlane> cutPl;

                cutPl->SetOrigin(arCutPos);
                cutPl->SetNormal(arDir);
                plCol->AddItem(cutPl);
                vtkNew<vtkPlaneCollection> plColNeg;
                vtkNew<vtkPlane> cutPlNeg;
                cutPlNeg->SetOrigin(arCutPos);
                cutPlNeg->SetNormal(negArDir);
                plColNeg->AddItem(cutPlNeg);
                auto negSide = lVtkHelper::clipClosedKeepCellLabels(curPoly, plColNeg);
                vecClipedPolys.push_back(negSide);
                // cut from bottom up -> the remainder upper part is the new part to cut
                curPoly = lVtkHelper::clipClosedKeepCellLabels(poly, plCol);
                // std::cout << "#pln: " << array2Str3d(arCutPos) << " dir: " << array2Str3d(arDir) << std::endl;
                /*if (negSide->GetNumberOfCells() == 0 || curPoly->GetNumberOfCells() == 0)
                {
                    std::cout << "pln: " << array2Str3d(arCutPos) << " dir: " << array2Str3d(arDir) << std::endl;
                    double newbnds[6];
                    curPoly->GetBounds(newbnds);
                    std::cout << "bnds: " << array2Str6d(bnds) << " newbdns: " << array2Str6d(newbnds) << std::endl;
                    lVtkHelper::writeVtkToVtk(curPoly, "/home/louis/CAD/debugGA/posSide.vtk");
                    lVtkHelper::writeVtkToVtk(negSide, "/home/louis/CAD/debugGA/negSide.vtk");
                }*/
            }
            vecClipedPolys.push_back(curPoly);
            return vecClipedPolys;
        }

        // create two planes (one per collection) to cut the mesh through the center of mass along the longes of the 3 main axes (x,y,z)
        // std::get<0>(tuple) is the positive direction, std::get<1>(tuple) the negative
        std::tuple<vtkSmartPointer<vtkPlaneCollection>, vtkSmartPointer<vtkPlaneCollection>> getCuttingPlanesTroughCenterOfMass(vtkSmartPointer<vtkPolyData> originalMesh)
        {
            // cout << "in the löoop" << endl;
            auto curPoly = originalMesh;

            double bnds[6];
            curPoly->GetBounds(bnds);
            bounds2DiffIp(bnds);
            auto dirMaxBnd = std::distance(bnds, std::max_element(bnds, bnds + 3));

            // calculate the center of mass
            vtkNew<vtkCenterOfMass> centerOMFltr;
            centerOMFltr->SetInputData(curPoly);
            centerOMFltr->SetUseScalarsAsWeights(false);
            centerOMFltr->Update();
            double center[3];
            centerOMFltr->GetCenter(center);

            double normal[3];
            double invNormal[3];
            switch (dirMaxBnd)
            {
            case 0:
                normal[0] = 1;
                normal[1] = 0;
                normal[2] = 0;
                invNormal[0] = -1;
                invNormal[1] = 0;
                invNormal[2] = 0;
                break;

            case 1:
                normal[0] = 0;
                normal[1] = 1;
                normal[2] = 0;
                invNormal[0] = 0;
                invNormal[1] = -1;
                invNormal[2] = 0;
                break;

            case 2:
                normal[0] = 0;
                normal[1] = 0;
                normal[2] = 1;
                invNormal[0] = 0;
                invNormal[1] = 0;
                invNormal[2] = -1;
                break;

            default:
                // could should never come here
                assert(0);
                normal[0] = 1;
                normal[1] = 0;
                normal[2] = 0;
                invNormal[0] = -1;
                invNormal[1] = 0;
                invNormal[2] = 0;
                break;
            }
            /*cout << "maxIdx was: " << dirMaxBnd << endl;
            cout << "values bnds were (dif[0:2], orig[3:5]) " << array2Str<double, 6>(bnds) << endl;
            cout << "center is: " << array2Str<double, 3>(center) << endl;
            cout << "normal is: " << array2Str<double, 3>(normal) << endl;*/

            vtkNew<vtkPlane> cutPlane;
            cutPlane->SetOrigin(center);
            cutPlane->SetNormal(normal);
            vtkNew<vtkPlane> invCutPlane;
            invCutPlane->SetOrigin(center);
            invCutPlane->SetNormal(invNormal);

            vtkNew<vtkPlaneCollection> planeCol;
            planeCol->AddItem(cutPlane);
            vtkNew<vtkPlaneCollection> invPlaneCol;
            invPlaneCol->AddItem(invCutPlane);
            return std::make_tuple<vtkSmartPointer<vtkPlaneCollection>, vtkSmartPointer<vtkPlaneCollection>>(planeCol, invPlaneCol);
        }

        // TODO: cut the part smaller until it is printable and return a vector of new objects
        // TODO: this is as described in the paper BUT DOES NOT GUARANTEE PRINTABILITY
        //  example given:
        //      the printer has a large and flat build volume, resulting in a large diagonal.
        //      the part is a cube with side lengths of more than the printers height. The diagonal might be smaller than the printers diagonal.
        //      but there is no way that the part can be rotated to make it printable
        //      -> every orientation results in the cube being higher than the printers max print height
        /**
         * @brief cut this mesh unitl every part does fit in the printers build volume, the part gets oriented to the print dir before all the checks and cutting gets done, and transformed to default afterwards.
         *
         * @return std::vector<MeshWithAttributes> the printable parts
         */
        std::vector<MeshWithAttributes> cutUntilPrintable()
        {
            // how this method works:
            // transform, cut, inverse transform

            auto originalPrintDir = this->printDir;
            bool hasCh = true;
            if (!this->convexHull || this->convexHull->GetNumberOfCells() <= 0)
            {
                hasCh = false;
            }
            // VPD by plane cut through center of mass
            // direction = ? maybe the longest side?
            // center of mass of the part?
            // for part in parts vector: if not feasable split until feasable
            // one vec not feasable, one feasable. Do until no unfeasable

            // is printable --> NOOOOOOOOOOOOOOOOOOOOOO
            /*if (checkFeasibility(this) == 0) {
                std::vector<MeshWithAttributes> vecThisMeshwa;
                vecThisMeshwa.push_back(*this);
                return vecThisMeshwa;
            }*/

            auto tr = getPrintDirAsTransform();
            auto invTr = getPrintDirAsTransform();
            tr->Inverse();
            vtkNew<vtkTransformPolyDataFilter> transformFltr;
            transformFltr->SetTransform(tr);
            transformFltr->SetInputData(this->mesh);
            transformFltr->Update();
            vtkSmartPointer<vtkPolyData> transfPoly = transformFltr->GetOutput();
            vtkNew<vtkTransformPolyDataFilter> transformFltrCh;
            transformFltrCh->SetTransform(tr);
            transformFltrCh->SetInputData(this->convexHull);
            transformFltrCh->Update();
            vtkSmartPointer<vtkPolyData> transfCh = transformFltrCh->GetOutput();

            std::vector<vtkSmartPointer<vtkPolyData>> vecUnfeasPoly;
            std::vector<vtkSmartPointer<vtkPolyData>> vecFeasPoly;
            vecUnfeasPoly.push_back(transfPoly);

            std::vector<vtkSmartPointer<vtkPolyData>> vecChUnfeasPoly;
            std::vector<vtkSmartPointer<vtkPolyData>> vecChFeasPoly;
            if (hasCh)
            {
                vecChUnfeasPoly.push_back(transfCh);
            }

            // cut perpendicualr to the longest side until all are feasable
            // with the use of vtkMassProperties, a cut perp to the dir with most volume could be cut
            while (vecUnfeasPoly.size() > 0)
            {
                auto curPoly = vecUnfeasPoly.at(vecUnfeasPoly.size() - 1);
                vecUnfeasPoly.pop_back();

                // if feasable do nothing
                if (checkFeasibility(curPoly) == 0)
                {
                    vecFeasPoly.push_back(curPoly);
                    if (hasCh)
                    {
                        auto curChPoly = vecChUnfeasPoly.at(vecChUnfeasPoly.size() - 1);
                        vecChUnfeasPoly.pop_back();
                        vecChFeasPoly.push_back(curChPoly);
                    }
                    continue;
                }

                // the get planes and cutting are done in two steps so the same planes can be easily used for the convex hull
                auto tupCutPlaneCols = this->getCuttingPlanesTroughCenterOfMass(curPoly);
                auto planeCol = std::get<0>(tupCutPlaneCols);
                auto invPlaneCol = std::get<1>(tupCutPlaneCols);

                // here one could improve performance by doing the updateCellDataWithPointData only once all cuts are done
                // would need adjustments to convertCellData2PointsData, so that the interface points array does not get overwritten but updated

                // so the convex hull becomes a 0 size hull -> bad and hows that psbl=?
                auto posClip = lVtkHelper::clipClosedKeepCellLabels(curPoly, planeCol);
                auto negClip = lVtkHelper::clipClosedKeepCellLabels(curPoly, invPlaneCol);
                vtkSmartPointer<vtkPolyData> posChClip;
                vtkSmartPointer<vtkPolyData> negChClip;
                if (hasCh)
                {
                    auto curChPoly = vecChUnfeasPoly.at(vecChUnfeasPoly.size() - 1);
                    vecChUnfeasPoly.pop_back();
                    posChClip = lVtkHelper::clipClosedKeepCellLabels(curChPoly, planeCol);
                    negChClip = lVtkHelper::clipClosedKeepCellLabels(curChPoly, invPlaneCol);
                }

                // this ones kinda bad as it shouldn't happen ...
                if (posClip->GetNumberOfCells() > 0)
                {
                    if (checkFeasibility(posClip) == 0)
                    {
                        vecFeasPoly.push_back(posClip);
                        if (hasCh)
                            vecChFeasPoly.push_back(posChClip);
                    }
                    else
                    {
                        vecUnfeasPoly.push_back(posClip);
                        if (hasCh)
                            vecChUnfeasPoly.push_back(posChClip);
                    }
                }
                if (negClip->GetNumberOfCells() > 0)
                {
                    if (checkFeasibility(negClip) == 0)
                    {
                        vecFeasPoly.push_back(negClip);
                        if (hasCh)
                            vecChFeasPoly.push_back(negChClip);
                    }
                    else
                    {
                        vecUnfeasPoly.push_back(negClip);
                        if (hasCh)
                            vecChUnfeasPoly.push_back(negChClip);
                    }
                }
            }
            // if one would use the as described in the paper version one would need to
            // return vecFeasPoly and also only cuts at most once

            // make all parts printable (unlike described in the paper) -> needed to ensure printability
            // split all directions until they would fit in the print volume
            // could be made smarter by checking differrent orientations but would also change cost
            std::vector<vtkSmartPointer<vtkPolyData>> vecPrintablePoly;
            // the ch is not really necessarily printable but the name should symbolise that it relates to the poly vec
            std::vector<vtkSmartPointer<vtkPolyData>> vecPrintableCh;
            typedef std::vector<vtkSmartPointer<vtkPolyData>>::size_type vec_size_type;

            // how it should look:
            /** z-dir:
             * cutUntilPrintable
             * x-y-dir:
             * align longest axes
             * x-dir:
             * cutUntilPrintable
             * y-dir:
             * cutUntilPrintable
             **/

            // fast track for printable parts
            std::vector<vtkSmartPointer<vtkPolyData>> vecUnprintablePoly;
            std::vector<vtkSmartPointer<vtkPolyData>> vecUnprintableCh;
            double bnds[6];
            for (vec_size_type i = 0; i < vecFeasPoly.size(); i++)
            {
                auto poly = vecFeasPoly.at(i);
                assert(poly->GetNumberOfCells());
                auto printBbx = this->printerBuildVolume;
                vtkBoundingBox partBbx;
                poly->GetBounds(bnds);
                partBbx.SetBounds(bnds);
                // fast track
                if (bbxFitsInBbx(partBbx, printBbx))
                {
                    vecPrintablePoly.push_back(poly);
                    if (hasCh)
                        vecPrintableCh.push_back(vecChFeasPoly.at(i));
                }
                else
                {
                    vecUnprintablePoly.push_back(poly);
                    if (hasCh)
                        vecUnprintableCh.push_back(vecChFeasPoly.at(i));
                }
            }
            // all printable parts got the fasttrack treatment

            // the z-part
            std::vector<vtkSmartPointer<vtkPolyData>> vecZCutPoly;
            std::vector<vtkSmartPointer<vtkPolyData>> vecZCutCh;
            for (vec_size_type i = 0; i < vecUnprintablePoly.size(); i++)
            {
                auto poly = vecUnprintablePoly.at(i);
                auto printBbx = this->printerBuildVolume;
                assert(poly->GetNumberOfCells());
                vtkBoundingBox partBbx;
                poly->GetBounds(bnds);
                partBbx.SetBounds(bnds);
                auto vecSubParts = cutAlongAxisUntilItFits(poly, partBbx.GetLength(2), printBbx.GetLength(2) * upperTol, 2);
                vecZCutPoly.insert(vecZCutPoly.end(), vecSubParts.begin(), vecSubParts.end());
                /*for (auto part: vecZCutPoly) {
                    std::cout << "nclls: " << part->GetNumberOfCells() << std::endl;
                }*/
                if (hasCh)
                {
                    auto vecSubChs = cutAlongAxisUntilItFits(vecUnprintableCh.at(i), partBbx.GetLength(2), printBbx.GetLength(2) * upperTol, 2);
                    vecZCutCh.insert(vecZCutCh.end(), vecSubChs.begin(), vecSubChs.end());
                }
            }

            // the x-y-part
            for (vec_size_type i = 0; i < vecZCutPoly.size(); i++)
            {
                auto poly = vecZCutPoly.at(i);
                auto printBbx = this->printerBuildVolume;
                vtkBoundingBox partBbx;
                poly->GetBounds(bnds);
                partBbx.SetBounds(bnds);
                if (poly->GetNumberOfCells() == 0)
                    continue;

                // check if flip might help:
                double partLens[3];
                double printLens[4];
                partBbx.GetLengths(partLens);
                printBbx.GetLengths(printLens);
                // 0, 0 -> no rot. 1, 1 -> no rot. 1, 0 or 1, 0 rot -> if sum == 1 -> rot
                int partMaxDir = (partLens[0] > partLens[1]) ? 0 : 1;
                int printMaxDir = (printLens[0] > printLens[1]) ? 0 : 1;
                bool useRot = (partMaxDir + printMaxDir) == 1;
                // "rotate" the print bbx before cutting
                vtkBoundingBox alignedPrintBbx;
                if (useRot)
                {
                    printBbx.GetBounds(bnds);
                    alignedPrintBbx.SetBounds(bnds[2], bnds[3], bnds[0], bnds[1], bnds[4], bnds[5]);
                }
                else
                {
                    alignedPrintBbx = printBbx;
                }

                // the x-part (returns size 1 arr if no cut needed)
                std::vector<vtkSmartPointer<vtkPolyData>> vecXZCutCh;
                auto vecXZCutPoly = cutAlongAxisUntilItFits(poly, partBbx.GetLength(0), alignedPrintBbx.GetLength(0) * upperTol, 0);
                if (hasCh)
                    vecXZCutCh = cutAlongAxisUntilItFits(vecZCutCh.at(i), partBbx.GetLength(0), alignedPrintBbx.GetLength(0) * upperTol, 0);

                // the y-part
                for (vec_size_type j = 0; j < vecXZCutPoly.size(); j++)
                {
                    poly = vecXZCutPoly.at(j);
                    if (poly->GetNumberOfCells() == 0)
                        continue;
                    poly->GetBounds(bnds);
                    partBbx.SetBounds(bnds);
                    auto vecYXZCutPoly = cutAlongAxisUntilItFits(poly, partBbx.GetLength(1), alignedPrintBbx.GetLength(1) * upperTol, 1);
                    vecPrintablePoly.insert(vecPrintablePoly.end(), vecYXZCutPoly.begin(), vecYXZCutPoly.end());
                    if (hasCh)
                    {
                        auto vecYXZCutCh = cutAlongAxisUntilItFits(vecXZCutCh.at(j), partBbx.GetLength(1), alignedPrintBbx.GetLength(1) * upperTol, 1);
                        vecPrintableCh.insert(vecPrintableCh.end(), vecYXZCutCh.begin(), vecYXZCutCh.end());
                    }
                }
            }

            std::vector<MeshWithAttributes> vecUntransformed;
            std::cout << "thread: " << std::this_thread::get_id() << " feas poly: " << vecFeasPoly.size() << " printable poly: " << vecPrintablePoly.size() << endl;

            // create the untransformed output data (all previous cuts were made on the oriented part)
            for (size_t i = 0; i < vecPrintablePoly.size(); i++)
            {
                MeshWithAttributes newMwa;
                newMwa.copyWithoutMeshes(*this);
                newMwa.printDir = originalPrintDir;
                auto curMesh = vecPrintablePoly.at(i);

                vtkNew<vtkTransformPolyDataFilter> reTransfFltr;
                reTransfFltr->SetTransform(invTr);
                reTransfFltr->SetInputData(curMesh);
                reTransfFltr->Update();
                // it should be enough to assign newMwa.mesh = reTransFltr... but ... it works now? so dont risk it
                vtkNew<vtkPolyData> msh;
                msh->DeepCopy(reTransfFltr->GetOutput());
                newMwa.mesh = msh;

                if (hasCh)
                {
                    auto curCh = vecPrintableCh.at(i);
                    vtkNew<vtkTransformPolyDataFilter> reTransfChFltr;
                    reTransfChFltr->SetTransform(invTr);
                    reTransfChFltr->SetInputData(curCh);
                    reTransfChFltr->Update();
                    reTransfChFltr->Update();
                    newMwa.convexHull = reTransfChFltr->GetOutput();
                }
                // here one could think about recalculating the print dir
                vecUntransformed.push_back(newMwa);
            }

            return vecUntransformed;
        }

        /**
         * @brief calculate the normal of the biggest cell of the convex hull as print direction (as described in the paper)
         *
         */
        void calculatePrintDirFromCh()
        {
            auto greatCellId = this->getGreatestCHCellId();
            vtkNew<vtkGenericCell> chCell;
            // auto chCell = this->convexHull->GetCell(greatCellId);
            this->convexHull->GetCell(greatCellId, chCell);
            double p0[3];
            double p1[3];
            double p2[3];
            chCell->GetPoints()->GetPoint(0, p0);
            chCell->GetPoints()->GetPoint(1, p1);
            chCell->GetPoints()->GetPoint(2, p2);
            double norm[3];
            vtkTriangle::ComputeNormal(p0, p1, p2, norm);
            auto tempNormVec = vtkVector3d(norm);
            auto invTempNormVec = lVtkHelper::invertVector(tempNormVec);
            this->printDir = invTempNormVec;
        }

        // calculates the area on call for tri cells
        double calculateMeshCellArea(vtkIdType i)
        {
            auto pts = this->mesh->GetCell(i)->GetPoints();
            return lVtkHelper::computeTriArea(pts);
        }

        /**
         * @brief Get the Greatest Facet Id of the object. The mesh needs to be a trimesh
         *
         * @return
         */
        void computeMeshGreatestCellId()
        {
            lVtkHelper::getGreatestMeshCellAndSize(this->mesh, this->greatestCellId, this->greatestCellSize);
        }

        void computeCHGreatestCellId()
        {
            lVtkHelper::getGreatestMeshCellAndSize(this->convexHull, this->greatestCHCellId, this->greatestCHCellSize);
        }

        // use vtkCleanPolydata to clean the mesh and replace the this->mesh pointer with the new datas pointer
        void cleanMesh()
        {
            this->mesh = lVtkHelper::getCleanPolydata(this->mesh);
        }

        /**
         * @brief append a cell array called Distances to this mesh. The values of the array are the distances from the cell center in normal direction until intersection with any other cell of this mesh
         *
         */
        void appendDistanceFieldIp()
        {
            // can all be replaced by a getter that calculates it if necessary else returns it
            // prequisites
            auto poly = this->mesh;
            auto cellCenters = this->getMeshCellCenters();
            auto cellNorms = this->getMeshCellNormals();
            double bnds[3];
            this->mesh->GetBounds(bnds);
            auto gapDist = bounds2DiagLen(bnds);
            // cout << "gap is: " << gapDist << endl;
            //  auto gapDist = this->gapDistance;
            //  used for raycasting
            //  auto ot = meshwa->getCellLocator();
            auto ot = this->getCellLocator();

            // do magic
            auto nCells = poly->GetNumberOfCells();

            vtkNew<vtkDoubleArray> dists;
            dists->SetNumberOfComponents(1);
            dists->SetNumberOfTuples(nCells);
            dists->SetName("Distances");

            double startOffset = .000001;
            auto dist = 0;
            double startPt[3];
            double offStartp[3];
            double endPt[3];

            for (vtkIdType i = 0; i < nCells; i++)
            {
                cellCenters->GetTuple(i, startPt);
                cellNorms->GetTuple(i, offStartp);
                cellNorms->GetTuple(i, endPt);
                // add the normalvector scaled to gapDist to the start point
                add3Ip(mul3Ip(endPt, gapDist), startPt);
                add3Ip(mul3Ip(offStartp, startOffset), startPt);
                // cout << " endp: " << array2Str<double, 3>(&endPt[0]) << endl;
                //  don't care for these right now, but one needs them for the method
                //  could be used for bbx creation so that in future only checks within bbx's need to be done
                vtkNew<vtkPoints> isPts;
                vtkNew<vtkIdList> isClls;
                ot->IntersectWithLine(offStartp, endPt, 0.0001 /*tol*/, isPts, isClls);

                // check if a valid collission occured
                if (isPts->GetNumberOfPoints() > 0)
                {
                    if (isClls->GetId(0) == i)
                    {
                        // collision is same cell
                        if (isPts->GetNumberOfPoints() > 1)
                        {
                            assert(isClls->GetId(1) != i);
                            dist = eNorm3(sub3Ip(&startPt[0], isPts->GetPoint(1)));
                        }
                        else
                        {
                            // only intersection with self
                            dist = gapDist;
                        }
                    }
                    else
                    {
                        dist = eNorm3(sub3Ip(&startPt[0], isPts->GetPoint(0)));
                    }
                }
                else
                {
                    // no point found ->distance is in inf but use some smaller for visual
                    dist = gapDist;
                }
                dists->SetValue(i, dist);
            }
            this->mesh->GetCellData()->AddArray(dists);
        }

        /**
         * @brief compute the cells that create to small gapps according to the gap criteria. this->gapCells
         *
         */
        void computeGapTris()
        {

            // can all be replaced by a getter that calculates it if necessary else returns it
            // prequisites
            this->createLinksIfNeeded();
            auto poly = this->mesh;
            auto cellCenters = this->getMeshCellCenters();
            auto cellNorms = this->getMeshCellNormals();
            auto gapD = this->gapDistance;
            auto cellAreas = this->getMeshTriAreas();
            // used for raycasting
            auto ot = this->getCellLocator();
            // auto ot = this->getObbTree(true);

            // do magic
            double startOffset = .0000001;
            double startPt[3];
            double offStartp[3];
            double endPt[3];
            double nOff[3];
            double nDist[3];
            auto gapCellsLst = vtkSmartPointer<vtkIdList>::New();
            auto pppt = vtkSmartPointer<vtkPoints>::New();
            this->gapPtsAndCost.gapCosts.clear();

            auto nCells = poly->GetNumberOfCells();
            // unsafe
            //#pragma omp parallel for
            for (vtkIdType i = 0; i < nCells; i++)
            {
                cellCenters->GetTuple(i, startPt);
                cellNorms->GetTuple(i, nOff);
                cellNorms->GetTuple(i, nDist);

                // create offset startpoint, so it wont intersect itself
                vtkMath::MultiplyScalar(nOff, startOffset);
                vtkMath::Add(startPt, nOff, offStartp);

                // create endpoint of raytrace / cast
                vtkMath::MultiplyScalar(nDist, gapD);
                vtkMath::Add(startPt, nDist, endPt);
                // add the normalvector scaled to gapDist to the start point

                // don't care for these right now, but one needs them for the method
                // could be used for bbx creation so that in future only checks within bbx's need to be done
                vtkNew<vtkPoints> isPts;
                vtkNew<vtkIdList> isClls;
                // 1 inside, -1 outside, 0 no intersect
                // tol is the tolerance for vtkCell intersect with line
                // tool seems to be used as if dist² < tol² line is inside cell
                // which of the following IntersectWithLines to use depens on the used locator / might not always be obvious or matter
                // idk if one's faster or not but the one that returns the point seems fine
                // auto isFound = ot->IntersectWithLine(offStartp, endPt, 0.0001 /*tol*/, isPts, isClls);
                // auto start = std::chrono::high_resolution_clock::now();
                auto isFound = ot->IntersectWithLine(offStartp, endPt, 0.0001 /*tol*/, isPts, isClls);
                // auto stop = std::chrono::high_resolution_clock::now();
                // auto dur = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
                // std::cout << "finding cell took µs: " << dur.count() << std::endl;
                if (isFound != 0)
                {
                    // this should never happen as it theoretically can't (based on the calulation of the offset startP)
                    if (isFound < 0)
                        std::cout << "isF s smaller 0 " << std::endl;
                    double arDd[3];
                    double x[3]; // the coord of intersect
                    vtkMath::Subtract(startPt, x, arDd);
                    gapCellsLst->InsertNextId(i);
                    this->gapPtsAndCost.gapCosts.push_back(cellAreas->GetComponent(i, 0));
                }
            }
            this->gapCells = gapCellsLst;
            this->gapPtsAndCost.gapStartAndIntersectPts = pppt;
        }
    };
}

#endif // MESH_WITH_ATTRIBUTES_HPP