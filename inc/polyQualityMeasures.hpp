#ifndef POLY_QUALITY_MEASURES_HPP
#define POLY_QUALITY_MEASURES_HPP

#include <vtkVector.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <unordered_set>
#include <vtkTriangle.h>
#include <vtkCell.h>
#include <vtkPoints.h>
#include <vtkOBBTree.h>
#include "meshWithAttributes.hpp"

namespace polyQualityMeasures {

inline bool printResults = false;

// get the id of the greatest facet of the input poly mesh, returns the size of the biggest cell
double getGreatestFacetId(MeshWithAttributes* meshwa);


/* #region quality measures as described in Oh et al. Part Decomposition and Evaluation Based on standard Design Guidelines for Additive Manufacturability and Assemblability */

// this is the roughness of one component. The full fomula summs up the cost of all components
double roughness(MeshWithAttributes* meshwa);

double overhangArea(MeshWithAttributes* meshwa);

// kinda more work
double sharpness(MeshWithAttributes* meshwa);

// this ones intensive -> raycasting "af" if not done right (ones for the whole model and then one could check if a part contains a critical region (check on in convex hull basis))
// calculate the center of all facets, then cast a ray in normal direction and check if it intersects the mesh /
// component somewhere in range [0, thresholdGap]. 
// If a intersection is found, add the area of the tri (caster) to the cost.
// this one's not really prod ready as it wastes a lot of power for something that could easily be done much more efficient
double gap(MeshWithAttributes* meshwa);

// not sure
// needs volume of convex hull and part -> gets calculated in vhacd -> check if addable to the convex hull
// or just recalculate
/**
 * @brief calculates the concavity of the mesh with respect to its convexHull.
 *  The MeshWithAttributs needs to contain a convexHull vtkPolydata.
 *  Both meshes need to be tri meshes and should contain no cutting faces, holes, ....
 *  The formula is volConvexHull - volMesh
 * @param meshwa 
 * @return double 
 */
double concavity(MeshWithAttributes* meshwa);

// check if the bbx diag is smaller than the printer diag
double feasibility(MeshWithAttributes* meshwa);
double feasibility(vtkSmartPointer<vtkPolyData> poly, double maxPrintDiagonal, double costFeasability);

// needs connectivity information or a set of 2d polydata / curve of the cut
// for cutting using vtkClipClosedSurface SetScalarModeToLables needs to be used
// this way interface faces can be identified by checking the lables array for 1's
// according to the paper the interface area needs to be multiplied by -2 -> needs to be done after the summation of all parts if area.
double interfaceArea(MeshWithAttributes* meshwa);

// needs some sort of meshwa collection or vector and just returns the size -> keep inline
double quantity(MeshWithAttributes* meshwa);

/* #endregion */


}

#endif // POLY_QUALITY_MEASURES_HPP