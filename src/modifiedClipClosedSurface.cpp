/**
 * @file modifiedClipClosedSurface.cpp
 * @author your name (you@domain.com)
 * @brief added parentId, the code is mostly copied from the one vtk uses
 * @version 0.1
 * @date 2022-08-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "modifiedClipClosedSurface.hpp"

#include <vtkPointData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>

vtkStandardNewMacro(ModClipClosed);



//------------------------------------------------------------------------------
// A helper class to quickly locate an edge, given the endpoint ids.
// It uses an stl map rather than a table partitioning scheme, since
// we have no idea how many entries there will be when we start.  So
// the performance is approximately log(n).

class vtkCCSEdgeLocatorNode
{
public:
  vtkCCSEdgeLocatorNode()
    : ptId0(-1)
    , ptId1(-1)
    , edgeId(-1)
    , next(nullptr)
  {
  }

  void FreeList()
  {
    vtkCCSEdgeLocatorNode* ptr = this->next;
    while (ptr)
    {
      vtkCCSEdgeLocatorNode* tmp = ptr;
      ptr = ptr->next;
      tmp->next = nullptr;
      delete tmp;
    }
  }

  vtkIdType ptId0;
  vtkIdType ptId1;
  vtkIdType edgeId;
  vtkCCSEdgeLocatorNode* next;
};

class vtkCCSEdgeLocator
{
private:
  typedef std::map<vtkIdType, vtkCCSEdgeLocatorNode> MapType;
  MapType EdgeMap;

public:
  static vtkCCSEdgeLocator* New() { return new vtkCCSEdgeLocator; }

  void Delete()
  {
    this->Initialize();
    delete this;
  }

  // Description:
  // Initialize the locator.
  void Initialize();

  // Description:
  // If edge (i0, i1) is not in the list, then it will be added and
  // a pointer for storing the new edgeId will be returned.
  // If edge (i0, i1) is in the list, then edgeId will be set to the
  // stored value and a null pointer will be returned.
  vtkIdType* InsertUniqueEdge(vtkIdType i0, vtkIdType i1, vtkIdType& edgeId);
};

void vtkCCSEdgeLocator::Initialize()
{
  for (MapType::iterator i = this->EdgeMap.begin(); i != this->EdgeMap.end(); ++i)
  {
    i->second.FreeList();
  }
  this->EdgeMap.clear();
}

vtkIdType* vtkCCSEdgeLocator::InsertUniqueEdge(vtkIdType i0, vtkIdType i1, vtkIdType& edgeId)
{
  // Ensure consistent ordering of edge
  if (i1 < i0)
  {
    vtkIdType tmp = i0;
    i0 = i1;
    i1 = tmp;
  }

  // Generate a integer key, try to make it unique
#ifdef VTK_USE_64BIT_IDS
  vtkIdType key = ((i1 << 32) ^ i0);
#else
  vtkIdType key = ((i1 << 16) ^ i0);
#endif

  vtkCCSEdgeLocatorNode* node = &this->EdgeMap[key];

  if (node->ptId1 < 0)
  {
    // Didn't find key, so add a new edge entry
    node->ptId0 = i0;
    node->ptId1 = i1;
    return &node->edgeId;
  }

  // Search through the list for i0 and i1
  if (node->ptId0 == i0 && node->ptId1 == i1)
  {
    edgeId = node->edgeId;
    return nullptr;
  }

  int i = 1;
  while (node->next != nullptr)
  {
    i++;
    node = node->next;

    if (node->ptId0 == i0 && node->ptId1 == i1)
    {
      edgeId = node->edgeId;
      return nullptr;
    }
  }

  // No entry for i1, so make one and return
  node->next = new vtkCCSEdgeLocatorNode;
  node = node->next;
  node->ptId0 = i0;
  node->ptId1 = i1;
  node->edgeId = static_cast<vtkIdType>(this->EdgeMap.size() - 1);
  return &node->edgeId;
}


// vtkCxxSetObjectMacro(ModClipClosed, ClippingPlanes, vtkPlaneCollection);

//------------------------------------------------------------------------------

int ModClipClosed::RequestData(vtkInformation* vtkNotUsed(request),
                               vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  //------------------------------------------------------------------------------
  // Get the info objects
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // Get the input and output
  vtkPolyData* input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Create objects needed for temporary storage
  if (this->IdList == nullptr)
  {
    this->IdList = vtkIdList::New();
  }

  // Get the input points
  vtkPoints* inputPoints = input->GetPoints();
  vtkIdType numPts = 0;
  int inputPointsType = VTK_FLOAT;
  if (inputPoints)
  {
    numPts = inputPoints->GetNumberOfPoints();
    inputPointsType = inputPoints->GetDataType();
  }

  // Force points to double precision, copy the point attributes
  vtkPoints* points = vtkPoints::New();
  points->SetDataTypeToDouble();
  points->SetNumberOfPoints(numPts);

  vtkPointData* pointData = vtkPointData::New();
  vtkPointData* inPointData = nullptr;

  if (this->PassPointData)
  {
    inPointData = input->GetPointData();
    pointData->InterpolateAllocate(inPointData, numPts, 0);
  }

  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double point[3];
    inputPoints->GetPoint(ptId, point);
    points->SetPoint(ptId, point);
    // Point data is not copied from input
    if (inPointData)
    {
      pointData->CopyData(inPointData, ptId, ptId);
    }
  }

  // An edge locator to avoid point duplication while clipping
  vtkCCSEdgeLocator* edgeLocator = vtkCCSEdgeLocator::New();

  // A temporary polydata for the contour lines that are triangulated
  vtkPolyData* tmpContourData = vtkPolyData::New();

  // The cell scalars
  vtkUnsignedCharArray* lineScalars = nullptr;
  vtkUnsignedCharArray* polyScalars = nullptr;
  vtkUnsignedCharArray* inputScalars = nullptr;

  // For input scalars: the offsets to the various cell types
  vtkIdType firstLineScalar = 0;
  vtkIdType firstPolyScalar = 0;
  vtkIdType firstStripScalar = 0;

  // Make the colors to be used on the data.
  int numberOfScalarComponents = 1;
  unsigned char colors[3][3];

  if (this->ScalarMode == VTK_CCS_SCALAR_MODE_COLORS)
  {
    numberOfScalarComponents = 3;
    vtkClipClosedSurface::CreateColorValues(
      this->BaseColor, this->ClipColor, this->ActivePlaneColor, colors);
  }
  else if (this->ScalarMode == VTK_CCS_SCALAR_MODE_LABELS)
  {
    colors[0][0] = 0;
    colors[1][0] = 1;
    colors[2][0] = 2;
  }

  // This is set if we have to work with scalars.  The input scalars
  // will be copied if they are unsigned char with 3 components, otherwise
  // new scalars will be generated.
  if (this->ScalarMode)
  {
    // Make the scalars
    lineScalars = vtkUnsignedCharArray::New();
    lineScalars->SetNumberOfComponents(numberOfScalarComponents);

    vtkDataArray* tryInputScalars = input->GetCellData()->GetScalars();
    // Get input scalars if they are RGB color scalars
    if (tryInputScalars && tryInputScalars->IsA("vtkUnsignedCharArray") &&
      numberOfScalarComponents == 3 && tryInputScalars->GetNumberOfComponents() == 3)
    {
      inputScalars = static_cast<vtkUnsignedCharArray*>(input->GetCellData()->GetScalars());

      vtkIdType numVerts = 0;
      vtkIdType numLines = 0;
      vtkIdType numPolys = 0;
      vtkCellArray* tmpCellArray = nullptr;
      if ((tmpCellArray = input->GetVerts()))
      {
        numVerts = tmpCellArray->GetNumberOfCells();
      }
      if ((tmpCellArray = input->GetLines()))
      {
        numLines = tmpCellArray->GetNumberOfCells();
      }
      if ((tmpCellArray = input->GetPolys()))
      {
        numPolys = tmpCellArray->GetNumberOfCells();
      }
      firstLineScalar = numVerts;
      firstPolyScalar = numVerts + numLines;
      firstStripScalar = numVerts + numLines + numPolys;
    }
  }

  // Break the input lines into segments, generate scalars for lines
  vtkCellArray* lines = vtkCellArray::New();
  if (input->GetLines() && input->GetLines()->GetNumberOfCells() > 0)
  {
    vtkClipClosedSurface::BreakPolylines(
      input->GetLines(), lines, inputScalars, firstLineScalar, lineScalars, colors[0]);
  }

  // Copy the polygons, convert strips to triangles
  vtkCellArray* polys = nullptr;
  int polyMax = 3;
  if ((input->GetPolys() && input->GetPolys()->GetNumberOfCells() > 0) ||
    (input->GetStrips() && input->GetStrips()->GetNumberOfCells() > 0))
  {
    // If there are line scalars, then poly scalars are needed too
    if (lineScalars)
    {
      polyScalars = vtkUnsignedCharArray::New();
      polyScalars->SetNumberOfComponents(numberOfScalarComponents);
    }

    polys = vtkCellArray::New();
    vtkClipClosedSurface::CopyPolygons(
      input->GetPolys(), polys, inputScalars, firstPolyScalar, polyScalars, colors[0]);
    vtkClipClosedSurface::BreakTriangleStrips(
      input->GetStrips(), polys, inputScalars, firstStripScalar, polyScalars, colors[0]);

    // Check if the input has polys and quads or just triangles
    vtkIdType npts = 0;
    const vtkIdType* pts = nullptr;
    vtkCellArray* inPolys = input->GetPolys();
    inPolys->InitTraversal();
    while (inPolys->GetNextCell(npts, pts))
    {
      if (npts > polyMax)
      {
        polyMax = npts;
      }
    }
  }

  // Get the clipping planes
  vtkPlaneCollection* planes = this->ClippingPlanes;

  // Arrays for storing the clipped lines and polys.
  vtkCellArray* newLines = vtkCellArray::New();
  vtkCellArray* newPolys = nullptr;
  if (polys)
  {
    newPolys = vtkCellArray::New();
  }

  // The point scalars, needed for clipping (not for the output!)
  vtkDoubleArray* pointScalars = vtkDoubleArray::New();

  // The line scalars, for coloring the outline
  vtkCellData* inLineData = vtkCellData::New();
  inLineData->CopyScalarsOn();
  inLineData->SetScalars(lineScalars);
  if (lineScalars)
  {
    lineScalars->Delete();
    lineScalars = nullptr;
  }

  // The poly scalars, for coloring the faces
  vtkCellData* inPolyData = vtkCellData::New();
  inPolyData->CopyScalarsOn();
  inPolyData->SetScalars(polyScalars);
  if (polyScalars)
  {
    polyScalars->Delete();
    polyScalars = nullptr;
  }

  // Also create output attribute data
  vtkCellData* outLineData = vtkCellData::New();
  outLineData->CopyScalarsOn();

  vtkCellData* outPolyData = vtkCellData::New();
  outPolyData->CopyScalarsOn();

  // Go through the clipping planes and clip the input with each plane
  vtkCollectionSimpleIterator iter;
  int numPlanes = 0;
  if (planes)
  {
    planes->InitTraversal(iter);
    numPlanes = planes->GetNumberOfItems();
  }

  vtkPlane* plane = nullptr;
  for (int planeId = 0; planes && (plane = planes->GetNextPlane(iter)); planeId++)
  {
    this->UpdateProgress((planeId + 1.0) / (numPlanes + 1.0));
    if (this->GetAbortExecute())
    {
      break;
    }

    // Is this the last cut plane?  If so, generate triangles.
    int triangulate = 5;
    if (planeId == numPlanes - 1)
    {
      triangulate = polyMax;
    }

    // Is this the active plane?
    int active = (planeId == this->ActivePlaneId);

    // Convert the plane into an easy-to-evaluate function
    double pc[4];
    plane->GetNormal(pc);
    pc[3] = -vtkMath::Dot(pc, plane->GetOrigin());

    // Create the clip scalars by evaluating the plane at each point
    vtkIdType numPoints = points->GetNumberOfPoints();
    pointScalars->SetNumberOfValues(numPoints);
    for (vtkIdType pointId = 0; pointId < numPoints; pointId++)
    {
      double p[3];
      points->GetPoint(pointId, p);
      double val = p[0] * pc[0] + p[1] * pc[1] + p[2] * pc[2] + pc[3];
      pointScalars->SetValue(pointId, val);
    }

    // Prepare the output scalars
    outLineData->CopyAllocate(inLineData, 0, 0);
    outPolyData->CopyAllocate(inPolyData, 0, 0);

    // Reset the locator
    edgeLocator->Initialize();

    // Clip the lines
    this->ClipLines(
      points, pointScalars, pointData, edgeLocator, lines, newLines, inLineData, outLineData);

    // Clip the polys
    if (polys)
    {
      // Get the number of lines remaining after the clipping
      vtkIdType numClipLines = newLines->GetNumberOfCells();

      // Cut the polys to generate more lines
      this->ClipAndContourPolys(points, pointScalars, pointData, edgeLocator, triangulate, polys,
        newPolys, newLines, inPolyData, outPolyData, outLineData);

      // Add scalars for the newly-created contour lines
      vtkUnsignedCharArray* scalars =
        vtkArrayDownCast<vtkUnsignedCharArray>(outLineData->GetScalars());

      if (scalars)
      {
        // Set the color to the active color if plane is active
        unsigned char* color = colors[1 + active];
        unsigned char* activeColor = colors[2];

        vtkIdType numLines = newLines->GetNumberOfCells();
        for (vtkIdType lineId = numClipLines; lineId < numLines; lineId++)
        {
          unsigned char oldColor[3];
          scalars->GetTypedTuple(lineId, oldColor);
          if (numberOfScalarComponents != 3 || oldColor[0] != activeColor[0] ||
            oldColor[1] != activeColor[1] || oldColor[2] != activeColor[2])
          {
            scalars->SetTypedTuple(lineId, color);
          }
        }
      }

      // Generate new polys from the cut lines
      vtkIdType cellId = newPolys->GetNumberOfCells();
      vtkIdType numClipAndContourLines = newLines->GetNumberOfCells();

      // Create a polydata for the lines
      tmpContourData->SetPoints(points);
      tmpContourData->SetLines(newLines);
      tmpContourData->BuildCells();

      this->TriangulateContours(
        tmpContourData, numClipLines, numClipAndContourLines - numClipLines, newPolys, pc);

      // Add scalars for the newly-created polys
      scalars = vtkArrayDownCast<vtkUnsignedCharArray>(outPolyData->GetScalars());

      if (scalars)
      {
        unsigned char* color = colors[1 + active];

        vtkIdType numCells = newPolys->GetNumberOfCells();
        if (numCells > cellId)
        {
          // The insert allocates space up to numCells-1
          scalars->InsertTypedTuple(numCells - 1, color);
          for (; cellId < numCells; cellId++)
          {
            scalars->SetTypedTuple(cellId, color);
          }
        }
      }

      // Add scalars to any diagnostic lines that added by
      // TriangulateContours().  In usual operation, no lines are added.
      scalars = vtkArrayDownCast<vtkUnsignedCharArray>(outLineData->GetScalars());

      if (scalars)
      {
        unsigned char color[3] = { 0, 255, 255 };

        vtkIdType numCells = newLines->GetNumberOfCells();
        if (numCells > numClipAndContourLines)
        {
          // The insert allocates space up to numCells-1
          scalars->InsertTypedTuple(numCells - 1, color);
          for (vtkIdType lineCellId = numClipAndContourLines; lineCellId < numCells; lineCellId++)
          {
            scalars->SetTypedTuple(lineCellId, color);
          }
        }
      }
    }

    // Swap the lines, points, etcetera: old output becomes new input
    vtkCellArray* tmp1 = lines;
    lines = newLines;
    newLines = tmp1;
    newLines->Initialize();

    if (polys)
    {
      vtkCellArray* tmp2 = polys;
      polys = newPolys;
      newPolys = tmp2;
      newPolys->Initialize();
    }

    vtkCellData* tmp4 = inLineData;
    inLineData = outLineData;
    outLineData = tmp4;
    outLineData->Initialize();

    vtkCellData* tmp5 = inPolyData;
    inPolyData = outPolyData;
    outPolyData = tmp5;
    outPolyData->Initialize();
  }

  // Delete the locator
  edgeLocator->Delete();

  // Delete the contour data container
  tmpContourData->Delete();

  // Delete the clip scalars
  pointScalars->Delete();

  // Get the line scalars
  vtkUnsignedCharArray* scalars = vtkArrayDownCast<vtkUnsignedCharArray>(inLineData->GetScalars());

  if (this->GenerateOutline)
  {
    output->SetLines(lines);
  }
  else if (scalars)
  {
    // If not adding lines to output, clear the line scalars
    scalars->Initialize();
  }

  // if false -> the output has no polys at all not even the input ones
  if (this->GenerateFaces)
  {
    if (polys) {
      // append legacy information (-1 -> no legacy)
      vtkNew<vtkIdTypeArray> parentCellId;
      parentCellId->SetName("ParentCellId");
      parentCellId->SetNumberOfComponents(1);
      parentCellId->SetNumberOfTuples(polys->GetNumberOfCells());
      parentCellId->Fill(-1);
      for (auto rel: this->newCellRelation) {
        // rel: <cellId, newCellId>
        auto cellId = std::get<0>(rel);
        auto newCellId = std::get<1>(rel);
        parentCellId->SetComponent(newCellId, 0 , cellId);
      }
      output->SetPolys(polys);
      output->GetCellData()->AddArray(parentCellId);
    } else {
      output->SetPolys(polys);
    }

    if (polys && scalars)
    {
      vtkUnsignedCharArray* pScalars =
        vtkArrayDownCast<vtkUnsignedCharArray>(inPolyData->GetScalars());

      vtkIdType m = scalars->GetNumberOfTuples();
      vtkIdType n = pScalars->GetNumberOfTuples();

      if (n > 0)
      {
        unsigned char color[3];
        color[0] = color[1] = color[2] = 0;

        // This is just to expand the array
        scalars->InsertTypedTuple(n + m - 1, color);

        // Fill in the poly scalars
        for (vtkIdType i = 0; i < n; i++)
        {
          pScalars->GetTypedTuple(i, color);
          scalars->SetTypedTuple(i + m, color);
        }
      }
    }
  }

  lines->Delete();

  if (polys)
  {
    polys->Delete();
  }

  if (this->ScalarMode == VTK_CCS_SCALAR_MODE_COLORS)
  {
    scalars->SetName("Colors");
    output->GetCellData()->SetScalars(scalars);
  }
  else if (this->ScalarMode == VTK_CCS_SCALAR_MODE_LABELS)
  {
    // Don't use VTK_UNSIGNED_CHAR or they will look like color scalars
    vtkSignedCharArray* categories = vtkSignedCharArray::New();
    categories->DeepCopy(scalars);
    categories->SetName("Labels");
    output->GetCellData()->SetScalars(categories);
    categories->Delete();
  }
  else
  {
    output->GetCellData()->SetScalars(nullptr);
  }

  newLines->Delete();
  if (newPolys)
  {
    newPolys->Delete();
  }

  inLineData->Delete();
  outLineData->Delete();
  inPolyData->Delete();
  outPolyData->Delete();

  // Finally, store the points in the output
  vtkClipClosedSurface::SqueezeOutputPoints(output, points, pointData, inputPointsType);
  output->Squeeze();

  points->Delete();
  pointData->Delete();

  return 1;

//---------------
}

void ModClipClosed::ClipAndContourPolys(vtkPoints *points, vtkDoubleArray *pointScalars,
                                        vtkPointData *pointData, vtkCCSEdgeLocator *edgeLocator, int triangulate,
                                        vtkCellArray *inputCells, vtkCellArray *outputPolys, vtkCellArray *outputLines,
                                        vtkCellData *inCellData, vtkCellData *outPolyData, vtkCellData *outLineData)
{

  vtkIdList *idList = this->IdList;

  // How many sides for output polygons?
  int polyMax = VTK_INT_MAX;
  if (triangulate)
  {
    if (triangulate < 4)
    { // triangles only
      polyMax = 3;
    }
    else if (triangulate == 4)
    { // allow triangles and quads
      polyMax = 4;
    }
  }

  int triangulationFailure = 0;

  // Go through all cells and clip them.
  vtkIdType numCells = inputCells->GetNumberOfCells();

  inputCells->InitTraversal();
  for (vtkIdType cellId = 0; cellId < numCells; cellId++)
  {
    vtkIdType numPts = 0;
    const vtkIdType *pts = nullptr;
    inputCells->GetNextCell(numPts, pts);
    idList->Reset();

    vtkIdType i1 = pts[numPts - 1];
    double v1 = pointScalars->GetValue(i1);
    int c1 = (v1 > 0);

    // The ids for the current edge: init j0 to -1 if i1 will be clipped
    vtkIdType j0 = (c1 ? i1 : -1);
    vtkIdType j1 = 0;

    // To store the ids of the contour line
    vtkIdType linePts[2];
    linePts[0] = 0;
    linePts[1] = 0;

    for (vtkIdType i = 0; i < numPts; i++)
    {
      // Save previous point info
      vtkIdType i0 = i1;
      double v0 = v1;
      int c0 = c1;

      // Generate new point info
      i1 = pts[i];
      v1 = pointScalars->GetValue(i1);
      c1 = (v1 > 0);

      // If at least one edge end point wasn't clipped
      if ((c0 | c1))
      {
        // If only one end was clipped, interpolate new point
        if ((c0 ^ c1))
        {
          vtkClipClosedSurface::InterpolateEdge(
              points, pointData, edgeLocator, this->Tolerance, i0, i1, v0, v1, j1);

          if (j1 != j0)
          {
            idList->InsertNextId(j1);
            j0 = j1;
          }

          // Save as one end of the contour line
          linePts[c0] = j1;
        }

        if (c1)
        {
          j1 = i1;

          if (j1 != j0)
          {
            idList->InsertNextId(j1);
            j0 = j1;
          }
        }
      }
    }

    // Insert the clipped poly
    vtkIdType numPoints = idList->GetNumberOfIds();

    if (numPoints > polyMax)
    {
      vtkIdType newCellId = outputPolys->GetNumberOfCells();

      // Triangulate the poly and insert triangles into output.
      if (!this->TriangulatePolygon(idList, points, outputPolys))
      {
        triangulationFailure = 1;
      }

      // Copy the attribute data to the triangle cells
      vtkIdType nCells = outputPolys->GetNumberOfCells();
      for (; newCellId < nCells; newCellId++)
      {
        outPolyData->CopyData(inCellData, cellId, newCellId);
        this->newCellRelation.push_back(std::make_tuple(cellId, newCellId));
      }
    }
    else if (numPoints > 2)
    {
      // Insert the polygon without triangulating it
      vtkIdType newCellId = outputPolys->InsertNextCell(idList);
      outPolyData->CopyData(inCellData, cellId, newCellId);
      this->newCellRelation.push_back(std::make_tuple(cellId, newCellId));
    }

    // Insert the contour line if one was created
    if (linePts[0] != linePts[1])
    {
      vtkIdType newCellId = outputLines->InsertNextCell(2, linePts);
      outLineData->CopyData(inCellData, cellId, newCellId);
    }
  }

  if (triangulationFailure && this->TriangulationErrorDisplay)
  {
    vtkErrorMacro("Triangulation failed, output may not be watertight");
  }

  // Free up the idList memory
  idList->Initialize();
}


