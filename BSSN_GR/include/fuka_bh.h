#pragma once
#include <cmath>
#include <functional>
#include "kadath_adapted_bh.hpp"
#include "Configurator/config_binary.hpp"
#include "coord_fields.hpp"
#include "bco_utilities.hpp"
#include "exporter_utilities.hpp"
#include <iostream>
#include <memory>
//#include <fftw3.h>

//! Macro to declare a pointer data member with associated trivial accessors.
//! Borrowed from Kadath master branch - credit S. Alain
#define ptr_data_member(type,identifier,smart_ptr_type) \
protected:\
    std:: smart_ptr_type##_ptr<type> identifier;\
public:\
    std:: smart_ptr_type##_ptr<type> const & get_##identifier() const {return identifier;}\
    std:: smart_ptr_type##_ptr<type> & get_##identifier() {return identifier;}

//! Macro to declare internal variable with a read-only trivial accessor.
#define internal_variable(type,identifier) \
protected:\
    type identifier;\
public:\
    optimal_access_type<type> get_##identifier () const {return identifier;}

/**
 * Template alias to select the (usually) optimal access type, either to pass an argument
 * or to return an object. \c cutoff_factor allows to let bigger objects be passed or
 * returned by value instead of reference to const.
 */
template<typename T,std::size_t cutoff_factor=1> using optimal_access_type =
    typename std::conditional<(sizeof(T) <= cutoff_factor*sizeof(T*)),T,T const &>::type;

namespace fuka_id_bh_importer {
using namespace Kadath;
using namespace Kadath::FUKA_Config;

class FUKA_BH_XCTS {
  public:
  using config_t = kadath_config_boost<BCO_BH_INFO>;

  public:
  // Global declaration of 'quants'
  std::vector<std::reference_wrapper<const Kadath::Scalar>> quants;

  // Global pointers to raw initial data solution!
  ptr_data_member(Kadath::Space_adapted_bh, space, unique);
  ptr_data_member(Kadath::Scalar, conf, unique);
  ptr_data_member(Kadath::Scalar, lapse, unique);
  ptr_data_member(Kadath::Vector, shift, unique);
  ptr_data_member(Kadath::Tensor, A, unique);

  internal_variable(double, adm_inf);
  internal_variable(double, coord_radius);
  internal_variable(double, hor_lapse);
  internal_variable(config_t, bconfig);
  internal_variable(std::vector<double>, quant_vals);
  
  private:
  // To fill the excision region we need to interpolate
  // the solution leading up to the excision surface
  // before extrapolating inside.
  int interpolation_order = 8;
  // Offset the start of interpolation from the excision radius
  // Not recommended for offset != 0
  double interpolation_offset = 0.;
  // radial spacing of interpolation points
  double delta_r_rel = 0.3;

  public:
  FUKA_BH_XCTS(std::string filename) : bconfig(config_t(filename)) {
    load_from_file();
    load_evolution_variables();
  }

  private:
  void load_from_file() {
    auto dat_filename = bconfig.space_filename();
    FILE* ff1 = fopen (dat_filename.c_str(), "r") ;
    
    space.reset(new Kadath::Space_adapted_bh{ff1});
    conf.reset(new Kadath::Scalar{*space, ff1});
    lapse.reset(new Kadath::Scalar{*space, ff1});
    shift.reset(new Kadath::Vector{*space, ff1});
    fclose(ff1);
  }

  private:
  // Function to load the metric
  void load_evolution_variables(){
    int ndom = space->get_nbr_domains();
    
    for (int i = 0; i < export_utils::NUM_VQUANTS; ++i)
      quants.push_back(std::cref(*conf));
    
    // Adding conformal factor and lapse to quants array
    quants[export_utils::PSI] = std::cref(*conf);
    quants[export_utils::ALP] = std::cref(*lapse);

    quants[export_utils::BETX] = std::cref((*shift)(1));
    quants[export_utils::BETY] = std::cref((*shift)(2));
    quants[export_utils::BETZ] = std::cref((*shift)(3));

    Base_tensor base(shift->get_basis());

    // reconstruct traceless extrinsic curvature
    Metric_flat fmet(*space, base);
    System_of_eqs syst(*space);
    fmet.set_system(syst, "f");

    syst.add_cst("PI", M_PI);
    syst.add_cst("P", *conf);
    syst.add_cst("N", *lapse);
    syst.add_cst("bet", *shift);
    syst.add_def(
      "A_ij = (D_i bet_j + D_j bet_i - 2. / 3.* D^k bet_k * f_ij) /2. / N");
    syst.add_def(ndom - 1, "Madm = -dr(P) / 2 / PI");
    adm_inf  = 
        space->get_domain(ndom-1)->integ(syst.give_val_def("Madm")()(ndom-1) , OUTER_BC);

    //Get horizon radius.  Needs to be space.BH<1,2>+2 since we're reading
    //r at the INNER_BC
    // Tensor A(syst.give_val_def("A"));
    // export_utils::add_tensor_refs(
    //   quants, {
    //     export_utils::AXX, 
    //     export_utils::AXY, 
    //     export_utils::AXZ, 
    //     export_utils::AYY, 
    //     export_utils::AYZ, 
    //     export_utils::AZZ
    //   }, A
    // );
    A.reset(new Kadath::Tensor{syst.give_val_def("A")});
    // Tensor A(syst.give_val_def("A"));
    Index ind(*A);
    
    quants[export_utils::AXX] = std::cref((*A)(ind));
    ind.inc();
    quants[export_utils::AXY] = std::cref((*A)(ind));
    ind.inc();
    quants[export_utils::AXZ] = std::cref((*A)(ind));
    ind.inc();
    ind.inc();
    quants[export_utils::AYY] = std::cref((*A)(ind));
    ind.inc();
    quants[export_utils::AYZ] = std::cref((*A)(ind));
    ind.inc();
    ind.inc();
    ind.inc();
    quants[export_utils::AZZ] = std::cref((*A)(ind));
    coord_radius = bco_utils::get_radius(space->get_domain(1),EQUI);
  }

  private:
  void interpolate_solution(double x, double y, double z){
    // Interpolated values at coordinate location
    quant_vals.resize(export_utils::NUM_VQUANTS);
    double x_loc = 0;
    double x_shifted = x - x_loc;
    double y_shifted = y;
    double const r2yz = y_shifted * y_shifted + z * z;
    double radius_from_bh = std::sqrt(x_shifted * x_shifted + r2yz);
    double rbh = coord_radius;

    // lambda function for filling excised region
    auto interp_f = [&](auto& ah_r, auto& extrap_r, auto& bh_ori) {
      // Avoid division by "0"
      extrap_r = (extrap_r <= 1e-7) ? 1e-7 : extrap_r;

      double theta = std::acos(z / extrap_r);
      
      // atan2 is needed here
      double phi = std::atan2(y_shifted, (x_shifted - bh_ori)); 

      // Where the filling takes places
      export_utils::spherical_turduck(
        quants, quant_vals, interpolation_order, delta_r_rel, interpolation_offset, 
        rbh, extrap_r, theta, phi, 2, bh_ori
      );
    };

    // Determine if we are inside the excision surface or not
    if(radius_from_bh <= (1. + interpolation_offset) * rbh) {
      interp_f(rbh, radius_from_bh, x_loc);
    } else {
      // If we are outside the excision surface, we can simply
      // interpolate the spectral solution

      // Construct Kadath point
      int ndim = 3;
      Kadath::Point abs_coords(ndim);
      abs_coords.set(1) = x_shifted;
      abs_coords.set(2) = y_shifted;
      abs_coords.set(3) = z;

      // Interpolate
      for (int k = 0; k < export_utils::NUM_VQUANTS; ++k) {
        quant_vals[k] = quants[k].get().val_point(abs_coords);
      }
    }
  }
  void compute_bssn_vals(double* vars) const {
    auto const psi = quant_vals[export_utils::PSI];
    auto const psi2 = psi * psi;
    auto const psi4 = psi2 * psi2;

    double gd[3][3];
    gd[0][0] = psi4;
    gd[0][1] = 0.0;
    gd[0][2] = 0.0;
    gd[1][0] = 0.0;
    gd[1][1] = psi4;
    gd[1][2] = 0.0;
    gd[2][0] = 0.0;
    gd[2][1] = 0.0;
    gd[2][2] = psi4;

    double Kd[3][3];
    Kd[0][0] = quant_vals[export_utils::AXX] * psi4;
    Kd[0][1] = quant_vals[export_utils::AXY] * psi4;
    Kd[0][2] = quant_vals[export_utils::AXZ] * psi4;
    Kd[1][0] = Kd[0][1];
    Kd[1][1] = quant_vals[export_utils::AYY] * psi4;
    Kd[1][2] = quant_vals[export_utils::AYZ] * psi4;
    Kd[2][0] = Kd[0][2];
    Kd[2][1] = Kd[1][2];
    Kd[2][2] = quant_vals[export_utils::AZZ] * psi4;

    double gu[3][3], gtd[3][3], Atd[3][3];
    double detgd, idetgd, trK;
    double chi;
    double t1, t2, t4, t6, t7, t9, t10, t12, t16;

    #include "adm2bssn.h"

    vars[VAR::U_SYMGT0] = gtd[0][0];
    vars[VAR::U_SYMGT1] = gtd[0][1];
    vars[VAR::U_SYMGT2] = gtd[0][2];
    vars[VAR::U_SYMGT3] = gtd[1][1];
    vars[VAR::U_SYMGT4] = gtd[1][2];
    vars[VAR::U_SYMGT5] = gtd[2][2];

    vars[VAR::U_SYMAT0] = Atd[0][0];
    vars[VAR::U_SYMAT1] = Atd[0][1];
    vars[VAR::U_SYMAT2] = Atd[0][2];
    vars[VAR::U_SYMAT3] = Atd[1][1];
    vars[VAR::U_SYMAT4] = Atd[1][2];
    vars[VAR::U_SYMAT5] = Atd[2][2];

    vars[VAR::U_K] = trK;
    vars[VAR::U_CHI] = chi;

    vars[VAR::U_BETA0] = 0.0;
    vars[VAR::U_BETA1] = 0.0;
    vars[VAR::U_BETA2] = 0.0;

    vars[VAR::U_GT0] = 0.0;
    vars[VAR::U_GT1] = 0.0;
    vars[VAR::U_GT2] = 0.0;

    vars[VAR::U_B0] = 0.0;
    vars[VAR::U_B1] = 0.0;
    vars[VAR::U_B2] = 0.0;
  }
  public:
  void operator()(double _x, double _y, double _z, double* var) { 
    interpolate_solution(_x, _y, _z);
    compute_bssn_vals(var); 
  }
};

// extern FUKA_BH_XCTS fuka_bh_data;
}