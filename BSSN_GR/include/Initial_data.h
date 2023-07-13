#pragma once
#include "grDef.h"
#include "parameters.h"
#include <memory>
#include "fuka_bh.h"
namespace initial_data
{
    // struct Initial_data {
        // inline void operator()(double _x, double _y, double _z, double* var) {
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
    // };
    // auto import_initial_data = [&](auto& fx) {
    //     return 
    // };
} // namespace initial_data

