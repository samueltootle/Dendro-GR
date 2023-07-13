#pragma once
#include <cmath>
#include <functional>
#include "kadath_adapted_bh.hpp"
#include "Configurator/config_binary.hpp"
#include "coord_fields.hpp"
#include "bco_utilities.hpp"
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
  enum ID_VARS{ 
    PSI, ALP, 
    BX   , BY   , BZ, 
    DXPSI, DYPSI, DZPSI, 
    DXALP, DYALP, DZALP, 
    DXBX , DYBX, DZBX, 
    DXBY , DYBY, DZBY, 
    DXBZ , DYBZ, DZBZ, 
    NUM_QUANTS 
  };

  public:
  // Global declaration of 'quants'
  std::vector<std::reference_wrapper<const Kadath::Scalar>> quants;

  // Global pointers to scalar & tensor arrays
  ptr_data_member(Kadath::Space_adapted_bh, space, unique);
  ptr_data_member(Kadath::Scalar, conf, unique);
  ptr_data_member(Kadath::Scalar, lapse, unique);
  ptr_data_member(Kadath::Vector, shift, unique);
  ptr_data_member(Kadath::Vector, grad_conf, unique);
  ptr_data_member(Kadath::Vector, grad_lapse, unique);
  ptr_data_member(Kadath::Vector, grad_betx, unique);
  ptr_data_member(Kadath::Vector, grad_bety, unique);
  ptr_data_member(Kadath::Vector, grad_betz, unique);

  internal_variable(double, adm_inf);
  internal_variable(double, coord_radius);
  internal_variable(double, hor_lapse);
  internal_variable(config_t, bconfig);

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
    
    // Adding conformal factor and lapse to quants array
    quants.push_back(std::cref(*conf));
    quants.push_back(std::cref(*lapse));

    auto add_all_components = [&](auto field) {
      for(size_t i = 1; i <=3; ++i)
        quants.push_back(std::ref((field)(i)));
    };

    // Add all components of shift to quants array
    add_all_components(*shift);
    
    grad_conf.reset(new Vector(conf->grad()));
    grad_lapse.reset(new Vector(lapse->grad()));
    grad_betx.reset(new Vector((*shift)(1).grad()));
    grad_bety.reset(new Vector((*shift)(2).grad()));
    grad_betz.reset(new Vector((*shift)(3).grad()));
    
    // Add all gradients to quants array
    add_all_components(*grad_conf);
    add_all_components(*grad_lapse);
    add_all_components(*grad_betx);
    add_all_components(*grad_bety);
    add_all_components(*grad_betz);

    Base_tensor base(shift->get_basis());

    // reconstruct traceless extrinsic curvature
    Metric_flat fmet(*space, base);
    System_of_eqs syst(*space);
    fmet.set_system(syst, "f");

    syst.add_cst("PI", M_PI);
    syst.add_cst("P", *conf);
    syst.add_def(ndom - 1, "Madm = -dr(P) / 2 / PI");
    adm_inf  = 
          space->get_domain(ndom-1)->integ(syst.give_val_def("Madm")()(ndom-1) , OUTER_BC);

    //Get horizon radius.  Needs to be space.BH<1,2>+2 since we're reading
    //r at the INNER_BC
    coord_radius = bco_utils::get_radius(space->get_domain(1),EQUI);
  }

  public:
  // get value for specific VAR - see enum
  double get_val(double _x, double _y, double _z, const int VAR){

        // construct Kadath point (ndim=3)
        Kadath::Point abs_coords(3);
        abs_coords.set(1) = _x;
        abs_coords.set(2) = _y;
        abs_coords.set(3) = _z;

        return quants[VAR].get().val_point(abs_coords);
  }
  void operator()(double _x, double _y, double _z, double* var) {}
};

// extern FUKA_BH_XCTS fuka_bh_data;
}