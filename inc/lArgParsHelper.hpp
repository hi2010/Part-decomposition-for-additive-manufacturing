/**
 * @file lArgParsHelper.hpp
 * @author iain, 0x90 https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c (input parser), rest is new, also input parser is modified
 * @brief 
 * @version 0.1
 * @date 2022-04-24
 * 
 * 
 */

#ifndef L_ARG_PARS_HELPER_HPP
#define L_ARG_PARS_HELPER_HPP

#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <type_traits>

#include "VHACD.h"
#include "workflows.hpp"
#include "lFileHelper.hpp"
#include "meshWithAttributes.hpp"

/**
 * @brief cli parser -> input to tokens / vars
 * 
 */
class InputParser{
    public:
        InputParser (const int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }
        /// @author iain
        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }
        /// @author iain
        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }

        template <typename T>
        T parsType(std::string option, T defaultVal) {
            std::string strVal = getCmdOption(option);
            if (!strVal.empty()) {
                if (std::is_same<T, double>::value) return stod(strVal);
                else if (std::is_same<T, long>::value) return stol(strVal);
                else if (std::is_same<T, unsigned int>::value) return static_cast<unsigned int>(stoul(strVal));
                else if (std::is_same<T, int>::value) return stoi(strVal);
            }
            return defaultVal;
        }

        std::vector<std::string> getTokens(){
            return this->tokens;
        }

        std::string getTokenString() {
            std::string resStr = "";
            for (auto tok: this->tokens) {
                resStr.append(tok+"\t");
            }
            return resStr;
        }

    private:
        std::vector <std::string> tokens;
};

/**
 * @brief show available cli options
 * TODO: add license printout
 */
inline void parsPrintHelp() {
    std::cout << "usage: [options] inputFile.3mf outputFilePath\fileStemName \n" \
    << "file stem name is any name like myResult which will get result specific naming appended e.g: myResult -> myResultoriresult.3mf" \
    << "style: --option {value}\n" \
    << "decompose part options are: \n" \
    << "--n_population --n_generations \n" \
    << "--cost_roughness --cost_overhang --cost_sharpness --cost_gap\n" \
    << "--cost_concavity --cost_feasibility --cost_interface\n" \
    << "--cost_quantity\n" \
    << "--printer_dimensions x_valxy_valxz_val in mm (buildspace cubic)\n\n" \
    << "all files on input and output are tried to be repaired\n" \
    << "it is possible, that the ga finishes before max n of generations is reached\n\n"
    << "--licenses will print the licenses of used libs"
    << std::endl;
}

inline void parsPrintLicenses() {
    
}


/**
 * @brief does nothing because vhacd params are not parsed
 * 
 * @param inpPars 
 * @return VHACD::IVHACD::Parameters 
 */
inline VHACD::IVHACD::Parameters parseVhacdParams(const InputParser &inpPars) {
    // nothing to do here
    VHACD::IVHACD::Parameters params;
    return params;
}

/**
 * @brief parses the optimization weights
 * 
 * @param inpPars 
 * @return sCostWeights 
 */
inline sCostWeights parseCostWeights(InputParser inpPars) {
    struct sCostWeights cw;
    cw.costRoughness = inpPars.parsType("--cost_roughness", 1.);
    cw.costOverhangArea = inpPars.parsType("--cost_overhang", 1.);
    cw.costSharpness = inpPars.parsType("--cost_sharpness", 50.);//1.);
    cw.costGap = inpPars.parsType("--cost_gap", 1.);
    cw.costConcavity = inpPars.parsType("--cost-concavity", .0000001);//1.);
    cw.costFeasibility = inpPars.parsType("--cost_feasibility", 1.);
    cw.costInterfaceArea = inpPars.parsType("--cost_interface", .001);//1.);
    cw.costQuantity = inpPars.parsType("--cost_quantity", 500.);//1.);
    return cw;
}

/**
 * @brief parse the optimization settings that can be set
 * 
 * @param inpPars 
 * @return GA_Type 
 */
inline GA_Type parseGASettings(InputParser inpPars) {
    GA_Type gares;
    // some options are set in the optimizeGA function (not the ones set here)

    gares.population = inpPars.parsType("--n_population", 20u);
    //
    gares.generation_max = inpPars.parsType("--n_generations", 1000);
    //
    return gares;
}

/**
 * @brief parse the dimensions of the printer build volume
 * 
 * @param inpPars 
 * @param meshwa build volume gets set on this object
 */
inline void parsePrintDimensions(InputParser inpPars, polyQualityMeasures::MeshWithAttributes& meshwa) {
    std::string defaultPrintDim = "no dim given";
    std::string sPrintDim = defaultPrintDim;
    if (inpPars.cmdOptionExists("--printer_dimensions")){
        sPrintDim = inpPars.getCmdOption("--printer_dimensions");
    }
    if (sPrintDim == defaultPrintDim) {
        std::cerr << "no printer dimension have been given. (--printer_dimensions x_val,y_val,z_val) The execution will stop because for decomposition printer dimensions are essential" << std::endl;
        std::cerr << "input args were: " << inpPars.getTokenString() << std::endl;
        exit(-1);
    }
    // was "x"
    std::string sepChar = ",";
    auto xPos0 = sPrintDim.find(sepChar);
    auto xPos1 = sPrintDim.find(sepChar, xPos0 + 1);
    double xLen, yLen, zLen;

    bool invalidPrinterDims = false;

    if ((xPos0 != std::string::npos) && (xPos1 != std::string::npos)) {
        xLen = stod(sPrintDim.substr(0, xPos0));
        yLen = stod(sPrintDim.substr(xPos0+1, xPos1-1));
        zLen = stod(sPrintDim.substr(xPos1+1));
        if (xLen <= 0 || yLen <= 0 || zLen <= 0) {
            invalidPrinterDims = true;
        }
    } else {
        invalidPrinterDims = true;
    }
    if (invalidPrinterDims) {
        std::cerr << "invalid printer dimensions, exiting : " << sPrintDim << std::endl;
        exit(-1);
    }

    meshwa.setBuildVolume(xLen, yLen, zLen);
    double lens[3];
    meshwa.printerBuildVolume.GetLengths(lens);
    double diagLen = meshwa.printerBuildVolume.GetDiagonalLength();
    std::cout << "printer build volume: \n" << "xlen, ylen, zlen: " << lens[0] << ", " << lens[1] << ", " << lens[2] \
    << "\ndiagonal:" << diagLen << std::endl;
}

/**
 * @brief run the program with parser / run the cli
 * 
 * @param inpPars InputParser object filled with the cli content (argc, argv)
 * @param argc the programs argc (from cli)
 * @param argv the programs argv (from cli)
 */
inline void runParsMode(const InputParser &inpPars, int argc, char* argv[]) {
    polyQualityMeasures::MeshWithAttributes meshwa;

    parsePrintDimensions(inpPars, meshwa);
    sCostWeights cw = parseCostWeights(inpPars);
    GA_Type gaSettings = parseGASettings(inpPars);

    struct OptimizationSettings optSettings;
    optSettings.costWts = cw;
    optSettings.gaSettings = gaSettings;

    // check if in file ex.
    std::string inFilePat = argv[argc-2];
    std::cout << "in file path is: " << inFilePat << std::endl;
    // check validity
    FILE* pFile = fopen64(inFilePat.c_str(), "r");
    if (pFile == NULL) {
        //throw std::runtime_error("Could not open input file");
        std::cerr << "Could not open input file, exiting" << std::endl;
        exit(-1);
    }
    // TODO: check if output path is possible and it should be a folder
    std::string outFilePat = argv[argc-1];
    std::cout << "out file path is: " << outFilePat << std::endl;

    CFileHelper fileHelp;
    meshwa.mesh = vtkSmartPointer<vtkPolyData>::New();
    meshwa.mesh = fileHelp.load3mfMeshAsVtkAndTryRepair(inFilePat);
    std::cout << "nDefects: " << fileHelp.checkForDefects(meshwa.mesh) << std::endl;

    std::cout << "starting optimization:\n" \
    << "input file path:\n" << inFilePat \
    << "output files path:\n" << outFilePat << std::endl;

    MyMiddleCost bestCost;
    VHACD::IVHACD::Parameters resParams = optimizeGA(meshwa, optSettings, &bestCost);
    std::cout << "\n\noptimization finished: (concavity=maxconcavity)\n" \
    << "writing results to:\n" << outFilePat \
    << "\noptimized parameters are:\n" \
    << "alpha=" << resParams.m_alpha << " beta=" << resParams.m_beta \
    << " concavity=" << resParams.m_concavity << "\n" \
    << "costs were:\n" \
    << bestCost.toStringCustomSep() << std::endl;
    splitMeshwaAndSave2File(&meshwa, resParams, outFilePat, "result", "hullResult", &bestCost);


}

#endif // L_ARG_PARS_HELPER_HPP
