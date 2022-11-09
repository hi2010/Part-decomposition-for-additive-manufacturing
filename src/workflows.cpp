#include <functional>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <mutex>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>

#include <vtkBoundingBox.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkFillHolesFilter.h>
#include <vtkDataSetSurfaceFilter.h>

#include "workflows.hpp"
#include "meshWithAttributes.hpp"
#include "polyQualityMeasures.hpp"
#include "VHACD.h"
#include "lVtkHelper.hpp"
#include "openGA.hpp"
#include "mesh3mfWriter.hpp"
#include "lVhacdHelper.hpp"

// TODO: for debug (mkdir)
#include <sys/stat.h>
#include <sstream>

using namespace polyQualityMeasures;

// TODO: prbly not a very well idea.
inline std::mutex mtx_vhacd_compute;

// calculate the whole cost for one part
// functions get only calculated if the cost is != 0
double calculateCost(MeshWithAttributes *meshwa, sCostWeights costs)
{
    // assumption: meshwa is the mesh and the convex hull. Init neededn't be done, as this is a per cost thing (some fct need some, some other not -> should be made in class of meshwa)
    double totalCost = 0;
    if (costs.costRoughness)
        totalCost += roughness(meshwa);
    if (costs.costOverhangArea)
        totalCost += overhangArea(meshwa);
    if (costs.costSharpness)
        totalCost += sharpness(meshwa);
    if (costs.costGap)
        totalCost += gap(meshwa);
    if (costs.costConcavity)
        totalCost += concavity(meshwa);
    // this ones rubbish if the splitting is already done until printable
    if (costs.costFeasibility)
        totalCost += feasibility(meshwa);
    if (costs.costInterfaceArea)
        totalCost += interfaceArea(meshwa);
    if (costs.costQuantity)
        totalCost += quantity(meshwa);
    return totalCost;
}

std::vector<MeshWithAttributes> computeConvexHulls(VHACD::IVHACD::Parameters params, std::vector<lVtkHelper::VHACD_POINT_TYPE> &vPoints,
                                                  std::vector<lVtkHelper::VHACD_INDEX_TYPE> &vTriangles)
{
    // lock to avoid strange locking problems from the static mutex in vhacd
    mtx_vhacd_compute.lock();

    // vhacd decompose and copy the objects to the vecMeshwa
    VHACD::IVHACD *ivhacd = VHACD::CreateVHACD();
    ivhacd->Compute(&vPoints[0], (unsigned int)vPoints.size() / 3,
                    (const uint32_t *)&vTriangles[0], (unsigned int)vTriangles.size() / 3, params);

    // do work on ivhacd (copy), so it can be freed and used by others
    auto nChs = ivhacd->GetNConvexHulls();

    std::vector<MeshWithAttributes> vecMeshwa;
    for (size_t i = 0; i < nChs; i++)
    {
        VHACD::IVHACD::ConvexHull ch;
        ivhacd->GetConvexHull(i, ch);
        MeshWithAttributes mwa;
        bool chSuccess;
        mwa.convexHull = vtkSmartPointer<vtkPolyData>::New();
        mwa.convexHull = lVtkHelper::convertVHACDCH2Polydata(ch, chSuccess);
        // this shouldn't happen but did so if an empty ch exists just ignore it
        if (chSuccess == false)
            continue;
        // only add meshes that contain data
        if (mwa.convexHull->GetNumberOfCells() > 0)
        {
            vecMeshwa.push_back(mwa);
        }
    }
    ivhacd->Clean();
    ivhacd->Release();
    mtx_vhacd_compute.unlock();
    return vecMeshwa;
}

std::vector<MeshWithAttributes> splitPoly(MeshWithAttributes *meshwa, VHACD::IVHACD::Parameters params)
{
    // convert mesh to (flat) vectors of points / indexes
    vtkSmartPointer<vtkPolyData> poly = meshwa->mesh;
    std::vector<lVtkHelper::VHACD_POINT_TYPE> vPoints;
    std::vector<lVtkHelper::VHACD_INDEX_TYPE> vTriangles;
    lVtkHelper::convertVtkPolyDataToFlatVectors(poly, vPoints, vTriangles);
    
    std::vector<MeshWithAttributes> vecMeshwa = computeConvexHulls(params, vPoints, vTriangles);
    std::cout << "thread: " << std::this_thread::get_id() << " nChs: " << vecMeshwa.size() << std::endl;

    // split the mesh until it is printable and add it to the parts vector if it contains cells
    std::vector<MeshWithAttributes> outVec;
    for (auto mwa : vecMeshwa)
    {
        // copy the attribs
        mwa.copyWithoutMeshes(*meshwa);
        vtkNew<vtkPolyData> tempPoly;
        tempPoly->DeepCopy(poly);
        mwa.mesh = tempPoly;
        
        mwa.calculatePrintDirFromCh();
        mwa.cutMeshWithCh();

        auto vecFeasParts = mwa.cutUntilPrintable();
        for (auto part : vecFeasParts)
        {
            if (part.mesh->GetNumberOfCells() > 0)
            {
                outVec.push_back(part);
            }
        }
    }

    return outVec;
}

std::vector<MeshWithAttributes> splitMeshwa(MeshWithAttributes *meshwa, VHACD::IVHACD::Parameters params)
{
    vtkSmartPointer<vtkPolyData> poly = meshwa->mesh;
    std::vector<lVtkHelper::VHACD_POINT_TYPE> vPoints;
    std::vector<lVtkHelper::VHACD_INDEX_TYPE> vTriangles;
    lVtkHelper::convertVtkPolyDataToFlatVectors(poly, vPoints, vTriangles);
    auto vecTemp = splitPoly(meshwa, params);
    return vecTemp;
}

/* target: give some file path some costs and some other config then let it run until finished
 the gen opt may be in another file (as class) this class may create the initial meshwa*/

/* #region copied partly from opnega example_so1.cpp */

double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X)
{
    // finalize the cost
    return X.middle_costs.totalCost;
}

std::ofstream output_file;

void SO_report_generation(
    int generation_number,
    const EA::GenerationType<MySolution, MyMiddleCost> &last_generation,
    const MySolution &best_genes)
{
    std::cout
        << "Generation [" << generation_number << "], "
        << "Best=" << last_generation.best_total_cost << ", "
        << "Average=" << last_generation.average_cost << ", "
        << "Best genes=(" << best_genes.to_string() << ")"
        << ", "
        << "Exe_time=" << last_generation.exe_time
        << std::endl;
}

void printNCells(MeshWithAttributes meshwa)
{
    std::cout << "in lamb " << meshwa.mesh->GetNumberOfCells() << std::endl;
}

// as i was unable to get it to work with an object in the struct i will try with an object containing the data
class MyClass
{
public:
    // size type for indexing
    typedef std::vector<MeshWithAttributes>::size_type typeVecSize;
    // the optimization weights
    sCostWeights _costWeights;
    // the mesh to decompose
    MeshWithAttributes *_meshwa;
    // the upper limit for the v-hacd parameters
    double upperLim = 1;
    MySolution bestGenes;
    MyMiddleCost bestCost;
    // unused
    double minChVol = -1;
    double maxChs = -1;
    // approximately the current population num in the generation (gets set to 0 on report generation and ++ on eval)
    unsigned int curPopul = 0;

    /**
     * @brief Construct a new My Class object
     *
     * @param costWeights the costs to use for the optimization
     * @param meshwa the mesh to use for the optimization
     */
    MyClass(const sCostWeights &costWeights, MeshWithAttributes meshwa) : _meshwa()
    {
        this->_costWeights = costWeights;
        this->_meshwa = new MeshWithAttributes(meshwa);
    }

    /**
     * @brief Construct a new My Class object / copy constructor
     *
     * @param mc the object to copy
     */
    MyClass(const MyClass &mc)
    {
        this->_costWeights = mc._costWeights;
        auto mwa = *mc._meshwa;
        this->_meshwa = new MeshWithAttributes(mwa);
    }

    ~MyClass() { delete this->_meshwa; }

    /**
     * @brief initialize the first generation of species
     *
     * @param p the class with the resulting v-hacd parameters
     * @param rnd01 the random value generator
     */
    void init_genes(MySolution &p, const std::function<double(void)> &rnd01)
    {
        auto rndInRng = [=]()
        { return upperLim * rnd01(); };
        // TODO: find reasonable ranges
        p.alpha = rndInRng();
        p.beta = rndInRng();
        p.gamma = rndInRng();
        p.maxConcavity = rndInRng();
    }

    /**
     * @brief decompose the object and evaluate the result
     *
     * @param p the v-hacd parameters to use
     * @param c the resulting middle cost
     * @return true always
     * @return false never
     */
    bool eval_solution(const MySolution &p, MyMiddleCost &c)
    {
        VHACD::IVHACD::Parameters params;
        params.m_alpha = p.alpha;
        params.m_beta = p.beta;
        params.m_concavity = p.maxConcavity;
        if (this->minChVol > 0)
        {
            params.m_minVolumePerCH = this->minChVol;
        }
        auto costs = this->_costWeights;

        MeshWithAttributes mwa;
        mwa.copyMwa(*this->_meshwa);
        auto vecMwa = splitMeshwa(&mwa, params);

        // calculate the weighted costs if the cost value is not 0, then add to the middle cost (not summed up cost)
        for (typeVecSize i = 0; i < vecMwa.size(); i++)
        {
            auto meshwa = &vecMwa.at(i);
            if (costs.costRoughness)
                c.costRoughness += costs.costRoughness * roughness(meshwa);
            if (costs.costOverhangArea)
                c.costOverhangArea += costs.costOverhangArea * overhangArea(meshwa);
            if (costs.costSharpness)
                c.costSharpness += costs.costSharpness * sharpness(meshwa);
            if (costs.costGap)
                c.costGap += costs.costGap * gap(meshwa);
            if (costs.costConcavity)
                c.costConcavity += costs.costConcavity * concavity(meshwa);
            // this ones rubbish if the splitting is already done until printable
            if (costs.costFeasibility)
                c.costFeasibility += costs.costFeasibility * feasibility(meshwa);
            // interface area is positive so a minus is used
            if (costs.costInterfaceArea)
                c.costInterfaceArea -= costs.costInterfaceArea * interfaceArea(meshwa);
            if (costs.costQuantity)
                c.costQuantity += costs.costQuantity * quantity(meshwa);
            // std::cout << "parts costs were: " << c.toString() << std::endl;
        }
        this->curPopul++;
        std::cout << "thread: " << std::this_thread::get_id() << "[npopul:" + std::to_string(this->curPopul) + "]: "
                  << "totalparts costs were: " << c.toString() << std::endl;
        // calculate the total cost as sum of all middle costs
        c.calculateTotalCost();
        return true;
    }

    MySolution mutate(const MySolution &X_base, const std::function<double(void)> &rnd01, double shrink_scale)
    {
        MySolution X_new;
        const double mu = 0.2 * shrink_scale; // mutation radius
        // this is crazy caus it might lock forever
        auto createNewRand = [=]()
        { return mu * fabs(rnd01() - rnd01()); };
        auto outOfRange = [=](auto a)
        { return a < 0 || a > upperLim; };
        auto createNewRandInRange = [=]()
        {
            auto newRand = 0.;
            do
            {
                newRand = createNewRand();
            } while (outOfRange(newRand));
            return newRand;
        };
        X_new = X_base;
        X_new.alpha += createNewRandInRange();
        X_new.beta += createNewRandInRange();
        X_new.gamma += createNewRandInRange();
        X_new.maxConcavity += createNewRandInRange();

        return X_new;
    }

    MySolution crossover(const MySolution &X1, const MySolution &X2, const std::function<double(void)> &rnd01)
    {
        MySolution X_new;
        auto crossGene = [=](auto geneValX, auto geneValY)
        {
            auto r = rnd01();
            return r * geneValX + (1.0 - r) * geneValY;
        };
        X_new.alpha = crossGene(X1.alpha, X2.alpha);
        X_new.beta = crossGene(X1.beta, X2.beta);
        X_new.gamma = crossGene(X1.gamma, X2.gamma);
        X_new.maxConcavity = crossGene(X1.maxConcavity, X2.maxConcavity);

        return X_new;
    }

    void SO_report_generation(
        int generation_number,
        const EA::GenerationType<MySolution, MyMiddleCost> &last_generation,
        const MySolution &best_genes)
    {
        this->curPopul = 0;
        std::cout
            << "Generation [" << generation_number << "], "
            << "Best=" << last_generation.best_total_cost << ", "
            << "Average=" << last_generation.average_cost << ", "
            << "Best genes=(" << best_genes.to_string() << ")"
            << ", "
            << "Exe_time=" << last_generation.exe_time
            << std::endl;

        this->bestGenes = best_genes;
        this->bestCost = last_generation.chromosomes.at(last_generation.best_chromosome_index).middle_costs;
    }
};

VHACD::IVHACD::Parameters optimizeGA(MeshWithAttributes meshwa)
{
    struct sCostWeights costw;
    output_file.open("./bin/result_so1.txt");
    output_file << "step"
                << "\t"
                << "x_best"
                << "\t"
                << "y_best"
                << "\t"
                << "cost_avg"
                << "\t"
                << "cost_best"
                << "\n";

    EA::Chronometer timer;
    timer.tic();

    std::cout << "in ga: " << meshwa.mesh->GetNumberOfCells() << std::endl;
    auto lmbb = [meshwa](std::string a)
    { return printNCells(meshwa); };
    lmbb("kk");

    // to keep the poly data one could prbly pass the eval function as part of an obj that contains the poly

    // assign poly data here
    meshwa.setBuildVolume(vtkBoundingBox(0, 200, 0, 150, 0, 200));

    costw.costConcavity = .0000001;
    costw.costInterfaceArea = .001;
    // costw.costGap = 0;
    // costw.costQuantity = 5000;
    costw.costQuantity = 100;
    costw.costSharpness = 50;

    std::cout << "meshwa sharp corners: " << meshwa.sharpCornersAngle << endl;
    MyClass myObj(costw, meshwa);
    std::cout << "weights: " << myObj._costWeights.costConcavity << std::endl;
    std::cout << "my obj exists"
              << " sharp cost: " << myObj._meshwa->sharpCornersAngle << std::endl;

    // test a min chv vol optimization reducing thin parts
    double plens[3];
    meshwa.printerBuildVolume.GetLengths(plens);
    auto minLen = *std::min_element(plens, plens + 3);
    // auto minVol = std::pow(minLen, 3)*.5;
    auto minVol = std::pow(minLen, 3) * .3;
    myObj.minChVol = minVol;
    double paBnds[6];
    meshwa.mesh->GetBounds(paBnds);
    auto pBbx = vtkBoundingBox(paBnds);
    double paLens[3];
    pBbx.GetLengths(paLens);
    auto maxPaLen = *std::max_element(paLens, paLens + 3);
    myObj.maxChs = (maxPaLen / minLen) * 4;

    // set ga settings
    GA_Type ga_obj;
    ga_obj.problem_mode = EA::GA_MODE::SOGA;
    // vhacd (performance) contains static variables e.g. s_depth which most likely get shared between threads causing undefined behav if multi obj exist. so no use for now ...
    ga_obj.multi_threading = true;
    // ga_obj.multi_threading=false;
    ga_obj.idle_delay_us = 1; // switch between threads quickly
    ga_obj.verbose = false;
    ga_obj.population = 20;
    ga_obj.generation_max = 1000;
    ga_obj.calculate_SO_total_fitness = calculate_SO_total_fitness;
    auto initGenesFct = [&myObj](MySolution &p, const std::function<double(void)> &rnd01)
    { myObj.init_genes(p, rnd01); };
    ga_obj.init_genes = initGenesFct;
    auto evalSolutionFct = [&myObj](const MySolution &p, MyMiddleCost &c)
    { return myObj.eval_solution(p, c); };
    ga_obj.eval_solution = evalSolutionFct;
    auto mutateFct = [&myObj](const MySolution &X_base, const std::function<double(void)> &rnd01, double shrink_scale)
    { return myObj.mutate(X_base, rnd01, shrink_scale); };
    ga_obj.mutate = mutateFct;
    auto crossoverFct = [&myObj](const MySolution &X1, const MySolution &X2, const std::function<double(void)> &rnd01)
    { return myObj.crossover(X1, X2, rnd01); };
    ga_obj.crossover = crossoverFct;
    auto reportFct = [&myObj](int generation_number, const EA::GenerationType<MySolution, MyMiddleCost> &last_generation, const MySolution &best_genes)
    { myObj.SO_report_generation(generation_number, last_generation, best_genes); };
    ga_obj.SO_report_generation = reportFct;
    ga_obj.best_stall_max = 10;
    // 20 % elite, if that is 0 but there are more than one species use 1
    // this is needed to avoid segfs from the ga "lib"
    unsigned int eliteCount = ga_obj.population * .2;
    eliteCount = (eliteCount == 0 && ga_obj.population > 1) ? 1 : eliteCount;
    ga_obj.elite_count = eliteCount;
    ga_obj.crossover_fraction = 0.7;
    ga_obj.mutation_rate = 0.4;
    ga_obj.solve();

    VHACD::IVHACD::Parameters resParams;
    resParams.m_alpha = myObj.bestGenes.alpha;
    resParams.m_beta = myObj.bestGenes.beta;
    // not in performance
    // resParams.m_gamma = myObj.bestGenes.gamma;
    resParams.m_concavity = myObj.bestGenes.maxConcavity;
    // retrieve cost here possible

    std::cout << "The problem is optimized in " << timer.toc() << " seconds." << std::endl;

    output_file.close();
    return resParams;
}

VHACD::IVHACD::Parameters optimizeGA(polyQualityMeasures::MeshWithAttributes meshwa, OptimizationSettings optSettings, MyMiddleCost *ptrMidCost)
{
    struct sCostWeights costw;
    output_file.open("./bin/result_so1.txt");
    output_file << "step"
                << "\t"
                << "x_best"
                << "\t"
                << "y_best"
                << "\t"
                << "cost_avg"
                << "\t"
                << "cost_best"
                << "\n";

    EA::Chronometer timer;
    timer.tic();

    MyClass myObj(optSettings.costWts, meshwa);

    // uses the printers build volume to define the min convex hull volume (optional and added so that small parts get avoided)
    // test a min chv vol optimization reducing thin parts
    double plens[3];
    meshwa.printerBuildVolume.GetLengths(plens);
    auto minLen = *std::min_element(plens, plens + 3);
    // auto minVol = std::pow(minLen, 3)*.5;
    auto minVol = std::pow(minLen, 3) * .3;
    myObj.minChVol = minVol;
    double paBnds[6];
    meshwa.mesh->GetBounds(paBnds);
    auto pBbx = vtkBoundingBox(paBnds);
    double paLens[3];
    pBbx.GetLengths(paLens);
    auto maxPaLen = *std::max_element(paLens, paLens + 3);
    myObj.maxChs = (maxPaLen / minLen) * 4;

    // TODO: there might be a bug that ga does not create the right amount of population (used 100 but only counted to 70, could also be some other problem like getting results earlier or so or some feature like don't redo the elite)

    GA_Type ga_obj = optSettings.gaSettings;
    ga_obj.problem_mode = EA::GA_MODE::SOGA;
    // vhacd (performance) contains static variables e.g. s_depth which most likely get shared between threads causing undefined behav if multi obj exist. so no use for now ...
    ga_obj.multi_threading = true;
    // ga_obj.multi_threading=false;
    ga_obj.idle_delay_us = 1; // switch between threads quickly
    ga_obj.verbose = false;
    ga_obj.calculate_SO_total_fitness = calculate_SO_total_fitness;
    auto initGenesFct = [&myObj](MySolution &p, const std::function<double(void)> &rnd01)
    { myObj.init_genes(p, rnd01); };
    ga_obj.init_genes = initGenesFct;
    auto evalSolutionFct = [&myObj](const MySolution &p, MyMiddleCost &c)
    { return myObj.eval_solution(p, c); };
    ga_obj.eval_solution = evalSolutionFct;
    auto mutateFct = [&myObj](const MySolution &X_base, const std::function<double(void)> &rnd01, double shrink_scale)
    { return myObj.mutate(X_base, rnd01, shrink_scale); };
    ga_obj.mutate = mutateFct;
    auto crossoverFct = [&myObj](const MySolution &X1, const MySolution &X2, const std::function<double(void)> &rnd01)
    { return myObj.crossover(X1, X2, rnd01); };
    ga_obj.crossover = crossoverFct;
    auto reportFct = [&myObj](int generation_number, const EA::GenerationType<MySolution, MyMiddleCost> &last_generation, const MySolution &best_genes)
    { myObj.SO_report_generation(generation_number, last_generation, best_genes); };
    ga_obj.SO_report_generation = reportFct;
    ga_obj.best_stall_max = 10;
    // 20 % elite, if that is 0 but there are more than one species use 1
    // this is needed to avoid segfs from the ga "lib"
    unsigned int eliteCount = ga_obj.population * .2;
    eliteCount = (eliteCount == 0 && ga_obj.population > 1) ? 1 : eliteCount;
    ga_obj.elite_count = eliteCount;
    ga_obj.crossover_fraction = 0.7;
    ga_obj.mutation_rate = 0.4;
    ga_obj.solve();

    VHACD::IVHACD::Parameters resParams;
    resParams.m_alpha = myObj.bestGenes.alpha;
    resParams.m_beta = myObj.bestGenes.beta;
    // not in performance
    // resParams.m_gamma = myObj.bestGenes.gamma;
    resParams.m_concavity = myObj.bestGenes.maxConcavity;
    // retrieve cost here possible
    if (ptrMidCost != nullptr)
    {
        // copy the values to the input (optional) pointer
        *ptrMidCost = myObj.bestCost;
    }

    std::cout << "The problem is optimized in " << timer.toc() << " seconds." << std::endl;

    output_file.close();
    return resParams;
}

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
void splitMeshwaAndSave2File(MeshWithAttributes *meshwa, VHACD::IVHACD::Parameters params, const std::string &path, const std::string &resPrefix, const std::string &hullPrefix, MyMiddleCost *thisCost)
{
    std::vector<vtkSmartPointer<vtkTransform>> vecPartTransfs;
    typedef std::vector<vtkSmartPointer<vtkTransform>>::size_type transfIdx;

    auto parts = splitMeshwa(meshwa, params);

    std::vector<vtkSmartPointer<vtkPolyData>>::size_type i = 0;

    // store every subpart and hull as 3mf and vtk in original and transformed rotation
    for (; i < parts.size(); i++)
    {
        auto part = parts.at(i);

        // "assembly" / transformation storage -> affine matrix
        auto partTransf = vtkSmartPointer<vtkTransform>::New();
        partTransf->Identity();

        // object storage
        auto oriMsh = part.getMeshOriginedRotatedToPrintDir(partTransf);
        //      assembly
        partTransf->Inverse();
        vecPartTransfs.push_back(partTransf);

        // not transformed vtk
        auto flpat = path + resPrefix + "_" + std::to_string(i) + ".vtk";
        // not transformed hull vtk
        auto hullpat = path + hullPrefix + "_" + std::to_string(i) + ".vtk";
        // transformed vtk
        auto oriflpat = path + "ori" + resPrefix + "_" + std::to_string(i) + ".vtk";
        // transformed 3mf
        auto ori3mfflpat = path + "ori" + resPrefix + "_" + std::to_string(i) + ".3mf";

        // clean the output mesh (try to close holes and other errors)
        vtkNew<vtkFillHolesFilter> fillHolesFtlr;
        fillHolesFtlr->SetHoleSize(100);
        fillHolesFtlr->SetInputData(oriMsh);
        vtkNew<vtkDataSetSurfaceFilter> surfFltr;
        surfFltr->SetInputConnection(fillHolesFtlr->GetOutputPort());
        surfFltr->Update();
        vtkNew<vtkGeometryFilter> geomFltr;
        geomFltr->SetInputConnection(surfFltr->GetOutputPort());
        geomFltr->Update();
        vtkSmartPointer<vtkPolyData> cleanOriMsh = geomFltr->GetOutput();

        lVtkHelper::writeVtkToVtk(part.mesh, flpat.c_str());
        lVtkHelper::writeVtkToVtk(part.convexHull, hullpat.c_str());
        lVtkHelper::writeVtkToVtk(cleanOriMsh, oriflpat.c_str());
        mesh3mfWriter::writeVtkPolyDataAs3mfModel(cleanOriMsh, ori3mfflpat);
    }

    // create the file with the transformation info
    auto transfPat = path + "assembly_" + resPrefix + ".txt";
    std::string assemblyStr;
    for (transfIdx j = 0; j < vecPartTransfs.size(); j++)
    {
        assemblyStr += "part: " + std::to_string(j) + "\n" + "transformation:\n";
        auto transf = vecPartTransfs.at(j);
        // 4x4 matrix
        // to csv with sep: "," and linebreak: ";"
        assemblyStr += lVtkHelper::getStringFromvtk4x4MatrixValues(transf);
    }
    std::ofstream assemblyFile(transfPat);
    assemblyFile << assemblyStr;
    assemblyFile.close();
    std::cout << "wrote assembly file to: " << transfPat << std::endl;

    // store vhacd parameters
    auto vhacdParamsPat = path + "vhacdparams_" + resPrefix + ".json";
    std::ofstream vhacdParamsFile(vhacdParamsPat);
    vhacdParamsFile << getParametersAsJsonString(params);
    vhacdParamsFile.close();
    std::cout << "wrote vhacd params file to: " << vhacdParamsPat << std::endl;

    // write the costs to a file if they were given to this func
    if (thisCost != nullptr)
    {
        auto costPat = path + "cost_" + resPrefix + ".json";
        auto costJson = thisCost->toStringJson();
        std::ofstream costFile(costPat);
        costFile << costJson;
        costFile.close();
        std::cout << "wrote costs to file: " << costPat << std::endl;
    }
}

/* #endregion copied partly from openga... */