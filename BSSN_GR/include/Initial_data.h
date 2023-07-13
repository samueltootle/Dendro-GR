#pragma once
#include "grDef.h"
#include "parameters.h"
#include <memory>
#include "fuka_bh.h"
#include "json.hpp"
using json = nlohmann::json;
namespace initial_data
{
    inline static void import(double _x, double _y, double _z, double* var) {
        const double xx=GRIDX_TO_X(_x);
        const double yy=GRIDY_TO_Y(_y);
        const double zz=GRIDZ_TO_Z(_z);
        
        #ifdef FUKA_BH_ID
            fuka_id_bh_importer::FUKA_BH_XCTS fuka_bh_data("test.info");
            fuka_bh_data(xx,yy,zz,var);
        #else
            bssn::punctureData(xx,yy,zz,var);
        #endif
    }
    inline static void extract_parameters(json& parfile) {

    }
} // namespace initial_data

