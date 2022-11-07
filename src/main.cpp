#define _CRT_SECURE_NO_WARNINGS
#include <iostream>

#include "lib3mf_implicit.hpp"


//#define _CRTDBG_MAP_ALLOC

#ifdef _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#include <stdlib.h>
#endif // _CRTDBG_MAP_ALLOC

#include "VHACD.h"
#include "oclHelper.h"
#include "lArgParsHelper.hpp"
#include <vtkSmartPointer.h>
#include <licensesAsStr.hpp>

#include <stdio.h>

#ifdef CL_VERSION_1_1
bool InitOCL(const unsigned int oclPlatformID, const unsigned int oclDeviceID, OCLHelper& oclHelper, std::ostringstream& msg);
#endif // CL_VERSION_1_1

// they work inplace
void modSmrtPtr(vtkSmartPointer<vtkTransform> tr) {
    std::cout << "in tha rot" << std::endl;
    tr->RotateZ(90);
    std::cout << "after the rot" << std::endl;
}


int main(int argc, char* argv[])
{
    // TODO: for now allow no arg mode
    InputParser inpPars(argc, argv);
    if ((argc <= 1)  ||  inpPars.cmdOptionExists("-h")) {
        parsPrintHelp();
        //return 0;
    } else if (inpPars.cmdOptionExists("--licenses"))
    {
        std::cout << LICENSES_TEXT << std::endl;
    } else {
        runParsMode(inpPars, argc, argv);
        return 0;
    }    
    return 0;
}
