/**
 * @file modifiedClipClosedSurface.hpp
 * @author your name (you@domain.com)
 * @brief a modified version of clip closed surface, that adds inheritance information in form of parent cell to the after clip parts. ParentCellId contains the id of the cell before the cut if it existed or the id of the cell from which this cell originated, i think.
 * @version 0.1
 * @date 2022-08-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef MODIFIED_CLIP_CLOSED_SURFACE_HPP
#define MODIFIED_CLIP_CLOSED_SURFACE_HPP

#include <vtkClipClosedSurface.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkPlaneCollection.h>
#include <vtkUnsignedCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolygon.h>
#include <vtkIdList.h>
#include <vtkSignedCharArray.h>

#include <vtkIncrementalOctreePointLocator.h>

#include <algorithm>
#include <map>
#include <utility>
#include <vector>
#include <tuple>

class vtkCCSEdgeLocator;

// modified so that new polys get an cellData array ParentCellId that has the number of the parent cell if one exists and -1 if none exists
// parent in this case means that the cell is the cell or a part of the cell in the input data with the given id 
// value of array at óutput cell (i) = inputArrayCellNumber
class ModClipClosed: public vtkClipClosedSurface {

    public:
    static ModClipClosed* New();
    vtkTypeMacro(ModClipClosed, vtkClipClosedSurface);

    

    int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
      vtkInformationVector* outputVector);

    void ClipAndContourPolys(vtkPoints* points, vtkDoubleArray* pointScalars,
  vtkPointData* pointData, vtkCCSEdgeLocator* edgeLocator, int triangulate,
  vtkCellArray* inputCells, vtkCellArray* outputPolys, vtkCellArray* outputLines,
  vtkCellData* inCellData, vtkCellData* outPolyData, vtkCellData* outLineData);

  // cellId, newCellId
  std::vector<std::tuple<vtkIdType, vtkIdType>> newCellRelation;
};


#endif