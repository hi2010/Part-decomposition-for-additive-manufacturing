/**
 * @file lVhacdHelper.cpp
 * @author your name (you@domain.com)
 * @brief copied from kmammou's v-hacd
 * @version 0.1
 * @date 2022-08-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <chrono>

#include "lVhacdHelper.hpp"
#include "VHACD.h"
#include "oclHelper.h"

using namespace VHACD;
using namespace std;

MyCallback::MyCallback(void) {}
MyCallback::~MyCallback(){};
void MyCallback::Update(const double overallProgress, const double stageProgress, const double operationProgress,
                        const char *const stage, const char *const operation)
{
    if (m_printUpdates)
        cout << setfill(' ') << setw(3) << (int)(overallProgress + 0.5) << "% "
             << "[ " << stage << " " << setfill(' ') << setw(3) << (int)(stageProgress + 0.5) << "% ] "
             << operation << " " << setfill(' ') << setw(3) << (int)(operationProgress + 0.5) << "%" << endl;
};

MyLogger::MyLogger(void) {}
MyLogger::MyLogger(const string &fileName) { OpenFile(fileName); }
MyLogger::~MyLogger(){};
void MyLogger::Log(const char *const msg)
{
    if (m_file.is_open())
    {
        m_file << msg;
        m_file.flush();
    }
    else if (m_printLog)
    {
        cout << "VHACD_Logger: " << msg << endl;
    }
}

void MyLogger::OpenFile(const string &fileName)
{
    m_file.open(fileName.c_str());
}

void Usage(const Parameters &params)
{
    std::ostringstream msg;
    msg << "V-HACD V" << VHACD_VERSION_MAJOR << "." << VHACD_VERSION_MINOR << endl;
    msg << "Syntax: testVHACD [options] --input infile.obj --output outfile.wrl --log logfile.txt" << endl
        << endl;
    msg << "Options:" << endl;
    msg << "       --input                     Wavefront .obj input file name" << endl;
    msg << "       --output                    VRML 2.0 output file name" << endl;
    msg << "       --log                       Log file name" << endl;
    msg << "       --resolution                Maximum number of voxels generated during the voxelization stage (default=100,000, range=10,000-16,000,000)" << endl;
    msg << "       --maxhulls                  Maximum number of convex hulls to produce." << endl;
    msg << "       --concavity                 Maximum allowed concavity (default=0.0025, range=0.0-1.0)" << endl;
    msg << "       --planeDownsampling         Controls the granularity of the search for the \"best\" clipping plane (default=4, range=1-16)" << endl;
    msg << "       --convexhullDownsampling    Controls the precision of the convex-hull generation process during the clipping plane selection stage (default=4, range=1-16)" << endl;
    msg << "       --alpha                     Controls the bias toward clipping along symmetry planes (default=0.05, range=0.0-1.0)" << endl;
    msg << "       --beta                      Controls the bias toward clipping along revolution axes (default=0.05, range=0.0-1.0)" << endl;
    msg << "       --gamma                     Controls the maximum allowed concavity during the merge stage (default=0.00125, range=0.0-1.0)" << endl;
    msg << "       --delta                     Controls the bias toward maximaxing local concavity (default=0.05, range=0.0-1.0)" << endl;
    msg << "       --pca                       Enable/disable normalizing the mesh before applying the convex decomposition (default=0, range={0,1})" << endl;
    msg << "       --mode                      0: voxel-based approximate convex decomposition, 1: tetrahedron-based approximate convex decomposition (default=0, range={0,1})" << endl;
    msg << "       --maxNumVerticesPerCH       Controls the maximum number of triangles per convex-hull (default=64, range=4-1024)" << endl;
    msg << "       --minVolumePerCH            Controls the adaptive sampling of the generated convex-hulls (default=0.0001, range=0.0-0.01)" << endl;
    msg << "       --convexhullApproximation   Enable/disable approximation when computing convex-hulls (default=1, range={0,1})" << endl;
    msg << "       --oclAcceleration           Enable/disable OpenCL acceleration (default=0, range={0,1})" << endl;
    msg << "       --oclPlatformID             OpenCL platform id (default=0, range=0-# OCL platforms)" << endl;
    msg << "       --oclDeviceID               OpenCL device id (default=0, range=0-# OCL devices)" << endl;
    msg << "       --help                      Print usage" << endl
        << endl;
    msg << "Examples:" << endl;
    msg << "       testVHACD.exe --input bunny.obj --output bunny_acd.wrl --log log.txt" << endl
        << endl;
    cout << msg.str();
    if (params.m_paramsVHACD.m_logger)
    {
        params.m_paramsVHACD.m_logger->Log(msg.str().c_str());
    }
}

void ParseParameters(int argc, char *argv[], Parameters &params)
{
    for (int i = 1; i < argc; ++i)
    {
        if (!strcmp(argv[i], "--input"))
        {
            if (++i < argc)
                params.m_fileNameIn = argv[i];
        }
        else if (!strcmp(argv[i], "--output"))
        {
            if (++i < argc)
                params.m_fileNameOut = argv[i];
        }
        else if (!strcmp(argv[i], "--log"))
        {
            if (++i < argc)
                params.m_fileNameLog = argv[i];
        }
        else if (!strcmp(argv[i], "--resolution"))
        {
            if (++i < argc)
                params.m_paramsVHACD.m_resolution = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--concavity"))
        {
            if (++i < argc)
                params.m_paramsVHACD.m_concavity = atof(argv[i]);
        }
        else if (!strcmp(argv[i], "--planeDownsampling"))
        {
            if (++i < argc)
                params.m_paramsVHACD.m_planeDownsampling = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--convexhullDownsampling"))
        {
            if (++i < argc)
                params.m_paramsVHACD.m_convexhullDownsampling = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--alpha"))
        {
            if (++i < argc)
                params.m_paramsVHACD.m_alpha = atof(argv[i]);
        }
        else if (!strcmp(argv[i], "--beta"))
        {
            if (++i < argc)
                params.m_paramsVHACD.m_beta = atof(argv[i]);
        }
        else if (!strcmp(argv[i], "--maxhulls"))
        {
            if (++i < argc)
                params.m_paramsVHACD.m_maxConvexHulls = atoi(argv[i]);
        }
        /*else if (!strcmp(argv[i], "--pca"))
        {
            if (++i < argc)
                params.m_paramsVHACD.m_pca = atoi(argv[i]);
        }*/
        /*else if (!strcmp(argv[i], "--mode"))
        {
            if (++i < argc)
                params.m_paramsVHACD.m_mode = atoi(argv[i]);
        }*/
        else if (!strcmp(argv[i], "--maxNumVerticesPerCH"))
        {
            if (++i < argc)
                params.m_paramsVHACD.m_maxNumVerticesPerCH = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--minVolumePerCH"))
        {
            if (++i < argc)
                params.m_paramsVHACD.m_minVolumePerCH = atof(argv[i]);
        }
        else if (!strcmp(argv[i], "--convexhullApproximation"))
        {
            if (++i < argc)
                params.m_paramsVHACD.m_convexhullApproximation = atoi(argv[i]);
        }
        /*else if (!strcmp(argv[i], "--oclAcceleration"))
        {
            if (++i < argc)
                params.m_paramsVHACD.m_oclAcceleration = atoi(argv[i]);
        }*/
        else if (!strcmp(argv[i], "--oclPlatformID"))
        {
            if (++i < argc)
                params.m_oclPlatformID = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--oclDeviceID"))
        {
            if (++i < argc)
                params.m_oclDeviceID = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--help"))
        {
            params.m_run = false;
        }
    }
    params.m_paramsVHACD.m_resolution = (params.m_paramsVHACD.m_resolution < 64) ? 0 : params.m_paramsVHACD.m_resolution;
    params.m_paramsVHACD.m_planeDownsampling = (params.m_paramsVHACD.m_planeDownsampling < 1) ? 1 : params.m_paramsVHACD.m_planeDownsampling;
    params.m_paramsVHACD.m_convexhullDownsampling = (params.m_paramsVHACD.m_convexhullDownsampling < 1) ? 1 : params.m_paramsVHACD.m_convexhullDownsampling;
}

#ifdef CL_VERSION_1_1
bool InitOCL(const unsigned int oclPlatformID, const unsigned int oclDeviceID, OCLHelper &oclHelper, std::ostringstream &msg)
{

    bool res = true;
    vector<string> info;
    res = oclHelper.GetPlatformsInfo(info, "\t\t");
    if (!res)
        return res;

    const size_t numPlatforms = info.size();
    msg << "\t Number of OpenCL platforms: " << numPlatforms << endl;
    for (size_t i = 0; i < numPlatforms; ++i)
    {
        msg << "\t OpenCL platform [" << i << "]" << endl;
        msg << info[i];
    }
    msg << "\t Using OpenCL platform [" << oclPlatformID << "]" << endl;
    res = oclHelper.InitPlatform(oclPlatformID);
    if (!res)
        return res;

    info.clear();
    res = oclHelper.GetDevicesInfo(info, "\t\t");
    if (!res)
        return res;

    const size_t numDevices = info.size();
    msg << "\t Number of OpenCL devices: " << numDevices << endl;
    for (size_t i = 0; i < numDevices; ++i)
    {
        msg << "\t OpenCL device [" << i << "]" << endl;
        msg << info[i];
    }
    msg << "\t Using OpenCL device [" << oclDeviceID << "]" << endl;
    res = oclHelper.InitDevice(oclDeviceID);
    return res;
}
#endif // CL_VERSION_1_1

void GetFileExtension(const string &fileName, string &fileExtension)
{
    size_t lastDotPosition = fileName.find_last_of(".");
    if (lastDotPosition == string::npos)
    {
        fileExtension = "";
    }
    else
    {
        fileExtension = fileName.substr(lastDotPosition, fileName.size());
        transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::toupper);
    }
}

void ComputeRandomColor(Material &mat)
{
    mat.m_diffuseColor[0] = mat.m_diffuseColor[1] = mat.m_diffuseColor[2] = 0.0f;
    while (mat.m_diffuseColor[0] == mat.m_diffuseColor[1] || mat.m_diffuseColor[2] == mat.m_diffuseColor[1] || mat.m_diffuseColor[2] == mat.m_diffuseColor[0])
    {
        mat.m_diffuseColor[0] = (rand() % 100) / 100.0f;
        mat.m_diffuseColor[1] = (rand() % 100) / 100.0f;
        mat.m_diffuseColor[2] = (rand() % 100) / 100.0f;
    }
}

bool LoadOFF(const string &fileName, vector<float> &points, vector<int> &triangles, IVHACD::IUserLogger &logger)
{
    FILE *fid = fopen(fileName.c_str(), "r");
    if (fid)
    {
        const string strOFF("OFF");
        char temp[1024];
        fscanf(fid, "%1023s", temp);
        if (string(temp) != strOFF)
        {
            logger.Log("Loading error: format not recognized \n");
            fclose(fid);
            return false;
        }
        else
        {
            int nv = 0;
            int nf = 0;
            int ne = 0;
            fscanf(fid, "%i", &nv);
            fscanf(fid, "%i", &nf);
            fscanf(fid, "%i", &ne);
            points.resize(nv * 3);
            triangles.resize(nf * 3);
            const int np = nv * 3;
            for (int p = 0; p < np; p++)
            {
                fscanf(fid, "%f", &(points[p]));
            }
            int s;
            for (int t = 0, r = 0; t < nf; ++t)
            {
                fscanf(fid, "%i", &s);
                if (s == 3)
                {
                    fscanf(fid, "%i", &(triangles[r++]));
                    fscanf(fid, "%i", &(triangles[r++]));
                    fscanf(fid, "%i", &(triangles[r++]));
                }
                else // Fix me: support only triangular meshes
                {
                    for (int h = 0; h < s; ++h)
                        fscanf(fid, "%i", &s);
                }
            }
            fclose(fid);
        }
    }
    else
    {
        logger.Log("Loading error: file not found \n");
        return false;
    }
    return true;
}

bool LoadOBJ(const string &fileName, vector<float> &points, vector<int> &triangles, IVHACD::IUserLogger &logger)
{
    const unsigned int BufferSize = 1024;
    FILE *fid = fopen(fileName.c_str(), "r");

    if (fid)
    {
        char buffer[BufferSize];
        int ip[4];
        float x[3];
        char *pch;
        char *str;
        while (!feof(fid))
        {
            if (!fgets(buffer, BufferSize, fid))
            {
                break;
            }
            else if (buffer[0] == 'v')
            {
                if (buffer[1] == ' ')
                {
                    str = buffer + 2;
                    for (int k = 0; k < 3; ++k)
                    {
                        pch = strtok(str, " ");
                        if (pch)
                            x[k] = (float)atof(pch);
                        else
                        {
                            return false;
                        }
                        str = NULL;
                    }
                    points.push_back(x[0]);
                    points.push_back(x[1]);
                    points.push_back(x[2]);
                }
            }
            else if (buffer[0] == 'f')
            {

                pch = str = buffer + 2;
                int k = 0;
                while (pch)
                {
                    pch = strtok(str, " ");
                    if (pch && *pch != '\n')
                    {
                        ip[k++] = atoi(pch) - 1;
                    }
                    else
                    {
                        break;
                    }
                    str = NULL;
                }
                if (k == 3)
                {
                    triangles.push_back(ip[0]);
                    triangles.push_back(ip[1]);
                    triangles.push_back(ip[2]);
                }
                else if (k == 4)
                {
                    triangles.push_back(ip[0]);
                    triangles.push_back(ip[1]);
                    triangles.push_back(ip[2]);

                    triangles.push_back(ip[0]);
                    triangles.push_back(ip[2]);
                    triangles.push_back(ip[3]);
                }
            }
        }
        fclose(fid);
    }
    else
    {
        logger.Log("File not found\n");
        return false;
    }
    return true;
}

bool SaveOFF(const string &fileName, const float *const &points, const int *const &triangles, const unsigned int &nPoints,
             const unsigned int &nTriangles, IVHACD::IUserLogger &logger)
{
    ofstream fout(fileName.c_str());
    if (fout.is_open())
    {
        size_t nV = nPoints * 3;
        size_t nT = nTriangles * 3;
        fout << "OFF" << std::endl;
        fout << nPoints << " " << nTriangles << " " << 0 << std::endl;
        for (size_t v = 0; v < nV; v += 3)
        {
            fout << points[v + 0] << " "
                 << points[v + 1] << " "
                 << points[v + 2] << std::endl;
        }
        for (size_t f = 0; f < nT; f += 3)
        {
            fout << "3 " << triangles[f + 0] << " "
                 << triangles[f + 1] << " "
                 << triangles[f + 2] << std::endl;
        }
        fout.close();
        return true;
    }
    else
    {
        logger.Log("Can't open file\n");
        return false;
    }
}

bool SaveVRML2(ofstream &fout, const double *const &points, const int *const &triangles,
               const unsigned int &nPoints, const unsigned int &nTriangles,
               const Material &material, IVHACD::IUserLogger &logger)
{
    if (fout.is_open())
    {
        fout.setf(std::ios::fixed, std::ios::floatfield);
        fout.setf(std::ios::showpoint);
        fout.precision(6);
        size_t nV = nPoints * 3;
        size_t nT = nTriangles * 3;
        fout << "#VRML V2.0 utf8" << std::endl;
        fout << "" << std::endl;
        fout << "# Vertices: " << nPoints << std::endl;
        fout << "# Triangles: " << nTriangles << std::endl;
        fout << "" << std::endl;
        fout << "Group {" << std::endl;
        fout << "    children [" << std::endl;
        fout << "        Shape {" << std::endl;
        fout << "            appearance Appearance {" << std::endl;
        fout << "                material Material {" << std::endl;
        fout << "                    diffuseColor " << material.m_diffuseColor[0] << " "
             << material.m_diffuseColor[1] << " "
             << material.m_diffuseColor[2] << std::endl;
        fout << "                    ambientIntensity " << material.m_ambientIntensity << std::endl;
        fout << "                    specularColor " << material.m_specularColor[0] << " "
             << material.m_specularColor[1] << " "
             << material.m_specularColor[2] << std::endl;
        fout << "                    emissiveColor " << material.m_emissiveColor[0] << " "
             << material.m_emissiveColor[1] << " "
             << material.m_emissiveColor[2] << std::endl;
        fout << "                    shininess " << material.m_shininess << std::endl;
        fout << "                    transparency " << material.m_transparency << std::endl;
        fout << "                }" << std::endl;
        fout << "            }" << std::endl;
        fout << "            geometry IndexedFaceSet {" << std::endl;
        fout << "                ccw TRUE" << std::endl;
        fout << "                solid TRUE" << std::endl;
        fout << "                convex TRUE" << std::endl;
        if (nV > 0)
        {
            fout << "                coord DEF co Coordinate {" << std::endl;
            fout << "                    point [" << std::endl;
            for (size_t v = 0; v < nV; v += 3)
            {
                fout << "                        " << points[v + 0] << " "
                     << points[v + 1] << " "
                     << points[v + 2] << "," << std::endl;
            }
            fout << "                    ]" << std::endl;
            fout << "                }" << std::endl;
        }
        if (nT > 0)
        {
            fout << "                coordIndex [ " << std::endl;
            for (size_t f = 0; f < nT; f += 3)
            {
                fout << "                        " << triangles[f + 0] << ", "
                     << triangles[f + 1] << ", "
                     << triangles[f + 2] << ", -1," << std::endl;
            }
            fout << "                ]" << std::endl;
        }
        fout << "            }" << std::endl;
        fout << "        }" << std::endl;
        fout << "    ]" << std::endl;
        fout << "}" << std::endl;
        return true;
    }
    else
    {
        logger.Log("Can't open file\n");
        return false;
    }
}

bool SaveOBJ(ofstream &fout, const double *const &points, const int *const &triangles, const unsigned int &nPoints,
             const unsigned int &nTriangles, const Material &material, IVHACD::IUserLogger &logger, int convexPart, int vertexOffset)
{
    if (fout.is_open())
    {

        fout.setf(std::ios::fixed, std::ios::floatfield);
        fout.setf(std::ios::showpoint);
        fout.precision(6);
        size_t nV = nPoints * 3;
        size_t nT = nTriangles * 3;

        fout << "o convex_" << convexPart << std::endl;

        if (nV > 0)
        {
            for (size_t v = 0; v < nV; v += 3)
            {
                fout << "v " << points[v + 0] << " " << points[v + 1] << " " << points[v + 2] << std::endl;
            }
        }
        if (nT > 0)
        {
            for (size_t f = 0; f < nT; f += 3)
            {
                fout << "f "
                     << triangles[f + 0] + vertexOffset << " "
                     << triangles[f + 1] + vertexOffset << " "
                     << triangles[f + 2] + vertexOffset << " " << std::endl;
            }
        }
        return true;
    }
    else
    {
        logger.Log("Can't open file\n");
        return false;
    }
}

void ComputeVhacdAndTimeIt(std::vector<lVtkHelper::VHACD_POINT_TYPE>& vPoints,
                           std::vector<lVtkHelper::VHACD_INDEX_TYPE>& vTriangles, VHACD::IVHACD::Parameters& params)
{
    // run V-HACD
    // After function call
    auto start = chrono::high_resolution_clock::now();
    cout << "computing V-HACD" << endl;
    IVHACD *interfaceVHACD = CreateVHACD();
    // bool res = interfaceVHACD->Compute(&points[0], (unsigned int)points.size() / 3,
    //     (const uint32_t *)&triangles[0], (unsigned int)triangles.size() / 3, params.m_paramsVHACD);
    bool res = interfaceVHACD->Compute(&vPoints[0], (unsigned int)vPoints.size() / 3,
                                       (const uint32_t *)&vTriangles[0], (unsigned int)vTriangles.size() / 3, params);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "V-HACD done: " << res << " took time [s]: " << duration.count() / pow(10, 6) << endl;
    interfaceVHACD->Clean();
    interfaceVHACD->Release();
}

std::string getParametersAsJsonString(VHACD::IVHACD::Parameters params) {
    if (VHACD_VERSION_MAJOR != 3 && VHACD_VERSION_MINOR != 1) {
        std::cerr << "getParametersAsJsonString is not defined for this version of VHACD, this might result in wrong or not working str converison" \
            << " used verion is: " << VHACD_VERSION_MAJOR << " " << VHACD_VERSION_MINOR << std::endl;
    }
    // the fist std::string() is c++ because it needs a string in the begining to concatenate to ... (((std::string + "nextStr") + "nextStr") ...)
    std::string resStr = std::string("{") \
        + "\n" + "\"VHACD_VERSION_MAJOR\": " + std::to_string(VHACD_VERSION_MAJOR) \
        + ",\n" + "\"VHACD_VERSION_MINOR\": " + std::to_string(VHACD_VERSION_MINOR) \
        + ",\n" + "\"alpha\": " + std::to_string(params.m_alpha) \
        + ",\n" + "\"beta\": " + std::to_string(params.m_beta) \
        + ",\n" + "\"concavity\": " + std::to_string(params.m_concavity) \
        + ",\n" + "\"minVolumePerCh\": " + std::to_string(params.m_minVolumePerCH) \
        + ",\n" + "\"resolution\": " + std::to_string(params.m_resolution) \
        + ",\n" + "\"maxNumVerticesPerCH\": " + std::to_string(params.m_maxNumVerticesPerCH) \
        + ",\n" + "\"planeDownsampling\": " + std::to_string(params.m_planeDownsampling) \
        + ",\n" + "\"convexHullDownsampling\": " + std::to_string(params.m_convexhullDownsampling) \
        + ",\n" + "\"convexHullApproximation\": " + std::to_string(params.m_convexhullApproximation) \
        + ",\n" + "\"maxConvexHulls\": " + std::to_string(params.m_maxConvexHulls) \
        + "\n}";
    return resStr;
}