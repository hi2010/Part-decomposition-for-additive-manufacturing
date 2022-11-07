#include <iostream>
#include <vector>

#include "l3mfHelper.hpp"

using namespace Lib3MF;

namespace l3mfHelper
{

    void printReaderWarnings(PReader reader)
    {
        // print all warnings that the reader created on fileread
        for (Lib3MF_uint32 iWarning = 0; iWarning < reader->GetWarningCount(); iWarning++)
        {
            Lib3MF_uint32 nErrorCode;
            std::string sWarningMessage = reader->GetWarning(iWarning, nErrorCode);
            std::cout << "Encountered warning #" << nErrorCode << " : " << sWarningMessage << std::endl;
        }
    }

    PModel getPModelFrom3mfFile(std::string sFilename, bool bPrintReaderWarnings/*=true*/)
    {
        // segmentation fault on file not accessible or whatever the lib3mf defines
        // other fault if library can't be found
        PWrapper wrapper = CWrapper::loadLibrary();
        // from here on the library loading worked
        PModel model = wrapper->CreateModel();
        PReader reader = model->QueryReader("3mf");
        reader->ReadFromFile(sFilename);
        if (bPrintReaderWarnings)
            printReaderWarnings(reader);
        return model;
    }

    void printMetadataGroup(PMetaDataGroup meta)
    {
        std::cout << "meta count: " << meta->GetMetaDataCount() << std::endl;
        for (size_t i = 0; i < meta->GetMetaDataCount(); i++)
        {
            auto cMeta = meta->GetMetaData(i);
            std::cout << std::to_string(i) << ": Meta: " << cMeta->GetName() << " type: " << cMeta->GetType()
                      << " must preserve: " << cMeta->GetMustPreserve()
                      << " key: " << cMeta->GetKey() << " value: " << cMeta->GetValue() << " nm spc " << cMeta->GetNameSpace() << std::endl;
        }
        std::cout << "meta finished" << std::endl;
    }

    void print3mfVersion(PWrapper wrapper)
    {
        Lib3MF_uint32 nMajor, nMinor, nMicro;
        wrapper->GetLibraryVersion(nMajor, nMinor, nMicro);
        std::cout << "lib3mf version = " << nMajor << "." << nMinor << "." << nMicro;
        std::string sReleaseInfo, sBuildInfo;
        if (wrapper->GetPrereleaseInformation(sReleaseInfo))
        {
            std::cout << "-" << sReleaseInfo;
        }
        if (wrapper->GetBuildInformation(sBuildInfo))
        {
            std::cout << "+" << sBuildInfo;
        }
        std::cout << std::endl;
    }

    std::vector<PMeshObject> getMeshObjsFrom3mfModel(PModel model)
    {
        std::vector<PMeshObject> vMshObjs;
        PObjectIterator objectIterator = model->GetObjects();
        while (objectIterator->MoveNext())
        {
            PObject object = objectIterator->GetCurrentObject();
            if (object->IsMeshObject())
            {
                auto mshObj = model->GetMeshObjectByID(object->GetResourceID());
                vMshObjs.push_back(mshObj);
            }
            else if (object->IsComponentsObject())
            {
                auto compObj = model->GetComponentsObjectByID(object->GetResourceID());
                std::cout << "object is component object #" << object->GetResourceID() << std::endl;
            }
            else
            {
                std::cout << "unknown / non std object #" << object->GetResourceID() << ": " << std::endl;
            }
        }
        return vMshObjs;
    }
}
