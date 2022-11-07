/**
 * @file vhacdHelper.hpp
 * @author your name (you@domain.com)
 * @brief This code is mostly copied from VHACD-Test from kmammou github: https://github.com/kmammou/v-hacd  (Day:11, Month:03, Year:2022)
 * @version 0.1
 * @date 2022-03-11
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef VHACD_HELPER_HPP
#define VHACD_HELPER_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "VHACD.h"
// TODO:
#include "lVtkHelper.hpp"

class MyCallback : public VHACD::IVHACD::IUserCallback {
public:
    MyCallback(void);
    ~MyCallback();
    void Update(const double overallProgress, const double stageProgress, const double operationProgress,
        const char* const stage, const char* const operation);
    bool m_printUpdates = true;
};

class MyLogger : public VHACD::IVHACD::IUserLogger {
public:
    MyLogger(void);
    MyLogger(const std::string& fileName);
    ~MyLogger();
    void Log(const char* const msg);
    void OpenFile(const std::string& fileName);

    // if true logs messages to std_out
    bool m_printLog = true;
private:
    std::ofstream m_file;
};

struct Material {

    float m_diffuseColor[3];
    float m_ambientIntensity;
    float m_specularColor[3];
    float m_emissiveColor[3];
    float m_shininess;
    float m_transparency;
    Material(void)
    {
        m_diffuseColor[0] = 0.5f;
        m_diffuseColor[1] = 0.5f;
        m_diffuseColor[2] = 0.5f;
        m_specularColor[0] = 0.5f;
        m_specularColor[1] = 0.5f;
        m_specularColor[2] = 0.5f;
        m_ambientIntensity = 0.4f;
        m_emissiveColor[0] = 0.0f;
        m_emissiveColor[1] = 0.0f;
        m_emissiveColor[2] = 0.0f;
        m_shininess = 0.4f;
        m_transparency = 0.5f;
    };
};

struct Parameters {
    unsigned int m_oclPlatformID;
    unsigned int m_oclDeviceID;
    std::string m_fileNameIn;
    std::string m_fileNameOut;
    std::string m_fileNameLog;
    bool m_run;
    VHACD::IVHACD::Parameters m_paramsVHACD;
    Parameters(void)
    {
        m_run = true;
        m_oclPlatformID = 0;
        m_oclDeviceID = 0;
        m_fileNameIn = "";
        m_fileNameOut = "output.wrl";
        m_fileNameLog = "log.txt";
    }
};

bool LoadOFF(const std::string& fileName, std::vector<float>& points, std::vector<int>& triangles, VHACD::IVHACD::IUserLogger& logger);
bool LoadOBJ(const std::string& fileName, std::vector<float>& points, std::vector<int>& triangles, VHACD::IVHACD::IUserLogger& logger);
bool SaveOFF(const std::string& fileName, const float* const& points, const int* const& triangles, const unsigned int& nPoints,
    const unsigned int& nTriangles, VHACD::IVHACD::IUserLogger& logger);
bool SaveVRML2(std::ofstream& fout, const double* const& points, const int* const& triangles, 
    const unsigned int& nPoints, const unsigned int& nTriangles, 
    const Material& material, VHACD::IVHACD::IUserLogger& logger);
bool SaveOBJ(std::ofstream& fout, const double* const& points, const int* const& triangles, const unsigned int& nPoints,
    const unsigned int& nTriangles, const Material& material, VHACD::IVHACD::IUserLogger& logger, int convexPart, int vertexOffset);
void GetFileExtension(const std::string& fileName, std::string& fileExtension);
void ComputeRandomColor(Material& mat);
void Usage(const Parameters& params);
void ParseParameters(int argc, char* argv[], Parameters& params);

void ComputeVhacdAndTimeIt(std::vector<lVtkHelper::VHACD_POINT_TYPE>& vPoints,
                           std::vector<lVtkHelper::VHACD_INDEX_TYPE>& vTriangles, VHACD::IVHACD::Parameters& params);

std::string getParametersAsJsonString(VHACD::IVHACD::Parameters params);

#endif  // VHACD_HELPER_HPP