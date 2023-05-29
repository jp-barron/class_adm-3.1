/** @file thermodynamics.h Documented includes for thermodynamics module */

#ifndef __THERMODYNAMICS__
#define __THERMODYNAMICS__

#include "background.h"
#include "evolver_ndf15.h"
#include "evolver_rkck.h"
#include "wrap_hyrec.h"
#include "wrap_recfast.h"
#include "injection.h"
//#include "lambertw.h"

/**
 * List of possible recombination algorithms.
 */

enum recombination_algorithm {
  recfast,
  hyrec
};

/**
 * List of possible reionization schemes.
 */

enum reionization_parametrization {
  reio_none,       /**< no reionization */
  reio_camb,       /**< reionization parameterized like in CAMB */
  reio_bins_tanh,  /**< binned reionization history with tanh inteprolation between bins */
  reio_half_tanh,  /**< half a tanh, instead of the full tanh */
  reio_many_tanh,  /**< similar to reio_camb but with more than one tanh */
  reio_inter       /**< linear interpolation between specified points */
};

/**
 * Is the input parameter the reionization redshift or optical depth?
 */

enum reionization_z_or_tau {
  reio_z,  /**< input = redshift */
  reio_tau /**< input = tau */
};

/**
 * Two useful smooth step functions, for smoothing transitions in recfast.
 */

#define f1(x) (-0.75*x*(x*x/3.-1.)+0.5)  /**< goes from 0 to 1 when x goes from -1 to 1 */
#define f2(x) (x*x*(0.5-x/3.)*6.)        /**< goes from 0 to 1 when x goes from  0 to 1 */

/**
 * All thermodynamics parameters and evolution that other modules need to know.
 *
 * Once initialized by thermodynamics_init(), contains all the necessary information on the thermodynamics, and in particular, a
 * table of thermodynamical quantities as a function of the redshift, used for interpolation in other modules.
 */

struct thermodynamics
{
  /** @name - input parameters initialized by user in input module (all other quantities are computed in this module, given these parameters
   *   and the content of the 'precision' and 'background' structures) */

  //@{

  double YHe;  /**< \f$ Y_{He} \f$: primordial helium mass fraction rho_He/(rho_H+rho_He),
                  close but not exactly equal to the density fraction 4*n_He/(n_H+4*n_He) */
  double bbn_alpha_sensitivity; /**< Related to variation of fundamental constants (sensitivity of YHe to alpha) */

  enum recombination_algorithm recombination; /**< recombination code */

  enum recfast_photoion_modes recfast_photoion_mode; /**< photo-ionization coefficient mode of the recfast algorithm */

  enum reionization_parametrization reio_parametrization; /**< reionization scheme */

  enum reionization_z_or_tau reio_z_or_tau; /**< is the input parameter the reionization redshift or optical depth? */

  double tau_reio; /**< if above set to tau, input value of reionization optical depth */

  double z_reio;   /**< if above set to z,   input value of reionization redshift */

  short compute_cb2_derivatives; /**< do we want to include in computation derivatives of baryon sound speed? */

  short compute_damping_scale; /**< do we want to compute the simplest analytic approximation to the photon damping (or diffusion) scale? */

  /** parameters for reio_camb */

  double reionization_width; /**< width of H reionization */

  double reionization_exponent; /**< shape of H reionization */

  double helium_fullreio_redshift; /**< redshift for of helium reionization */

  double helium_fullreio_width; /**< width of helium reionization */

  /** parameters for reio_bins_tanh */

  int binned_reio_num; /**< with how many bins do we want to describe reionization? */

  double * binned_reio_z; /**< central z value for each bin */

  double * binned_reio_xe; /**< imposed \f$ X_e(z)\f$ value at center of each bin */

  double binned_reio_step_sharpness; /**< sharpness of tanh() step interpolating between binned values */

  /** parameters for reio_many_tanh */

  int many_tanh_num; /**< with how many jumps do we want to describe reionization? */

  double * many_tanh_z; /**< central z value for each tanh jump */

  double * many_tanh_xe; /**< imposed \f$ X_e(z)\f$ value at the end of each jump (ie at later times)*/

  double many_tanh_width; /**< sharpness of tanh() steps */

    /** parameters for reio_inter */

  int reio_inter_num; /**< with how many jumps do we want to describe reionization? */

  double * reio_inter_z; /**< discrete z values */

  double * reio_inter_xe; /**< discrete \f$ X_e(z)\f$ values */

  /** parameters for energy injection */

  short has_exotic_injection; /**< true if some exotic mechanism
                                  injects energy and affects the
                                  evolution of ionization and/or
                                  temperature and/or other
                                  thermodynamics variables that are
                                  relevant for the calculation of CMB
                                  anisotropies (and spectral
                                  distorsions if requested). */

  struct injection in; /**< structure to store exotic energy injections and their energy deposition */

  double annihilation; /**< parameter describing CDM annihilation (f <sigma*v> / m_cdm, see e.g. 0905.0003) */

  short has_on_the_spot; /**< flag to specify if we want to use the on-the-spot approximation **/

  double decay; /**< parameter describing CDM decay (f/tau, see e.g. 1109.6322)*/

  double annihilation_variation; /**< if this parameter is non-zero,
				     the function F(z)=(f <sigma*v> /
				     m_cdm)(z) will be a parabola in
				     log-log scale between zmin and
				     zmax, with a curvature given by
				     annihlation_variation (must be
				     negative), and with a maximum in
				     zmax; it will be constant outside
				     this range */

  double annihilation_z; /**< if annihilation_variation is non-zero,
			     this is the value of z at which the
			     parameter annihilation is defined, i.e.
			     F(annihilation_z)=annihilation */

  double annihilation_zmax; /**< if annihilation_variation is non-zero,
				redshift above which annihilation rate
				is maximal */

  double annihilation_zmin; /**< if annihilation_variation is non-zero,
				redshift below which annihilation rate
				is constant */

  double annihilation_f_halo; /**< takes the contribution of DM annihilation in halos into account*/
  double annihilation_z_halo; /**< characteristic redshift for DM annihilation in halos*/

  double a_idm_dr;      /**< strength of the coupling between interacting dark matter and interacting dark radiation (idm-idr) */
  double b_idr;         /**< strength of the self coupling for interacting dark radiation (idr-idr) */
  double nindex_idm_dr; /**< temperature dependence of the interaction between dark matter and dark radiation */
  double m_idm_dr;      /**< dark matter mass for idm_dr */

  /** parameters for varying fundamental constants */

  short has_varconst; /**< presence of varying fundamental constants? */
    
  /* BEGIN #TWIN SECTOR */
  short has_twin;
  /* END TWIN SECTOR */

  /** START #TWIN SECTOR */
  double YHe_twin;
  double z_reco_twin;/** Approximate redshift below which H recombination takes place */
  double z_threshold_twin;
  /** END TWIN SECTOR */

  //@}

  /** @name - all indices for the vector of thermodynamical (=th) quantities stored in table */

  //@{

  int index_th_xe;            /**< ionization fraction \f$ x_e \f$ */
  int index_th_dkappa;        /**< Thomson scattering rate \f$ d \kappa / d \tau\f$ (units 1/Mpc) */
  int index_th_tau_d;         /**< Baryon drag optical depth */
  int index_th_ddkappa;       /**< scattering rate derivative \f$ d^2 \kappa / d \tau^2 \f$ */
  int index_th_dddkappa;      /**< scattering rate second derivative \f$ d^3 \kappa / d \tau^3 \f$ */
  int index_th_exp_m_kappa;   /**< \f$ exp^{-\kappa} \f$ */
  int index_th_g;             /**< visibility function \f$ g = (d \kappa / d \tau) * exp^{-\kappa} \f$ */
  int index_th_dg;            /**< visibility function derivative \f$ (d g / d \tau) \f$ */
  int index_th_ddg;           /**< visibility function second derivative \f$ (d^2 g / d \tau^2) \f$ */
  int index_th_dmu_idm_dr;    /**< scattering rate of idr with idm_dr (i.e. idr opacity to idm_dr scattering) (units 1/Mpc) */
  int index_th_ddmu_idm_dr;   /**< derivative of this scattering rate */
  int index_th_dddmu_idm_dr;  /**< second derivative of this scattering rate */
  int index_th_dmu_idr;       /**< idr self-interaction rate */
  int index_th_tau_idm_dr;    /**< optical depth of idm_dr (due to interactions with idr) */
  int index_th_tau_idr;       /**< optical depth of idr (due to self-interactions) */
  int index_th_g_idm_dr;      /**< visibility function of idm_idr */
  int index_th_cidm_dr2;      /**< interacting dark matter squared sound speed \f$ c_{dm}^2 \f$ */
  int index_th_Tidm_dr;       /**< temperature of DM interacting with DR \f$ T_{idm_dr} \f$ */
  int index_th_Tb;            /**< baryon temperature \f$ T_b \f$ */
  int index_th_dTb;           /**< derivative of baryon temperature */
  int index_th_wb;            /**< baryon equation of state parameter \f$ w_b = k_B T_b / \mu \f$ */
  int index_th_cb2;           /**< squared baryon adiabatic sound speed \f$ c_b^2 \f$ */
  int index_th_dcb2;          /**< derivative wrt conformal time of squared baryon sound speed \f$ d [c_b^2] / d \tau \f$ (only computed if some non-minimal tight-coupling schemes is requested) */
  int index_th_ddcb2;         /**< second derivative wrt conformal time of squared baryon sound speed  \f$ d^2 [c_b^2] / d \tau^2 \f$ (only computed if some non0-minimal tight-coupling schemes is requested) */
  int index_th_rate;          /**< maximum variation rate of \f$ exp^{-\kappa}\f$, g and \f$ (d g / d \tau) \f$, used for computing integration step in perturbation module */
  int index_th_r_d;           /**< simple analytic approximation to the photon comoving damping scale */

  int th_size;                /**< size of thermodynamics vector */

  /** START #TWIN SECTOR */
  /* TWIN: Thermodynamics indices and table */
  int index_th_xe_twin;         /**< ionization fraction \f$ x_e \f$ */
  int index_th_Tb_twin;         /**< baryon temperature \f$ T_b \f$ */
  int index_th_wb_twin;         /**< baryon equation of state parameter \f$ w_b \f$ */
  int index_th_cb2_twin;        /**< squared baryon adiabatic sound speed \f$ c_b^2 \f$ */
  int index_th_dkappa_twin; /**< Thomson scattering rate \f$ d \kappa / d \tau \f$ (units 1/Mpc) */
  int index_th_tau_d_twin;         /**< Baryon drag optical depth */
  int index_th_ddkappa_twin;       /**< scattering rate derivative \f$ d^2 \kappa / d \tau^2 \f$ */
  int index_th_dddkappa_twin;      /**< scattering rate second derivative \f$ d^3 \kappa / d \tau^3 \f$ */
  int index_th_exp_m_kappa_twin;   /**< \f$ exp^{-\kappa} \f$ */
  int index_th_g_twin;             /**< visibility function \f$ g = (d \kappa / d \tau) * exp^{-\kappa} \f$ */
  int index_th_dg_twin;            /**< visibility function derivative \f$ (d g / d \tau) \f$ */
  int index_th_ddg_twin;           /**< visibility function second derivative \f$ (d^2 g / d \tau^2) \f$ */
  int index_th_rate_twin;          /**< maximum variation rate of \f$ exp^{-\kappa}\f$, g and \f$ (d g / d \tau) \f$, used for computing integration step in perturbation module */
  /** END TWIN SECTOR */

  /** START #TWIN SECTOR */
  /*Added by Jared Barron- Twin recombination constants that used to be in the recombination structure */
  /*double YHe_twin;*/
  double Tnow_twin;
  double fHe_twin;
  double Nnow_twin;
  double CR_twin;
  double CB1_twin;
  double CB1_He1_twin;
  double CB1_He2_twin;
  double x_e_Lalpha;
  double CT_twin;    
   /** END TWIN SECTOR */

  //@}

  /** @name - thermodynamics interpolation tables */

  //@{

  int tt_size;                   /**< number of lines (redshift steps) in the tables */
  double * z_table;              /**< vector z_table[index_z] with values of redshift (vector of size tt_size) */
  double * tau_table;            /**< vector tau_table[index_tau] with values of conformal time (vector of size tt_size) */
  double * thermodynamics_table; /**< table thermodynamics_table[index_z*pth->tt_size+pba->index_th] with all other quantities (array of size th_size*tt_size) */

  //@}

  /** @name - table of their second derivatives, used for spline interpolation */

  //@{

  double * d2thermodynamics_dz2_table; /**< table d2thermodynamics_dz2_table[index_z*pth->tt_size+pba->index_th] with values of \f$ d^2 t_i / dz^2 \f$ (array of size th_size*tt_size) */

  //@}

  /** @name - characteristic quantities like redshift, conformal time and sound horizon at recombination */

  //@{

  double z_rec;   /**< z at which the visibility reaches its maximum (= recombination redshift) */
  double tau_rec; /**< conformal time at which the visibility reaches its maximum (= recombination time) */
  double rs_rec;  /**< comoving sound horizon at recombination */
  double ds_rec;  /**< physical sound horizon at recombination */
  double ra_rec;  /**< conformal angular diameter distance to recombination */
  double da_rec;  /**< physical angular diameter distance to recombination */
  double rd_rec;  /**< comoving photon damping scale at recombination */

  double z_star;  /**< redshift at which photon optical depth crosses one */
  double tau_star;/**< confirmal time at which photon optical depth crosses one */
  double rs_star; /**< comoving sound horizon at z_star */
  double ds_star; /**< physical sound horizon at z_star */
  double ra_star;  /**< conformal angular diameter distance to z_star */
  double da_star;  /**< physical angular diameter distance to z_star */
  double rd_star;  /**< comoving photon damping scale at z_star */

  double z_d;     /**< baryon drag redshift */
  double tau_d;   /**< baryon drag time */
  double ds_d;    /**< physical sound horizon at baryon drag */
  double rs_d;    /**< comoving sound horizon at baryon drag */

  double tau_cut; /**< at at which the visibility goes below a fixed fraction of the maximum visibility, used for an approximation in perturbation module */

  double angular_rescaling;      /**< [ratio ra_rec / (tau0-tau_rec)]: gives CMB rescaling in angular space relative to flat model (=1 for curvature K=0) */
  double tau_free_streaming;     /**< minimum value of tau at which free-streaming approximation can be switched on */
  double tau_idr_free_streaming; /**< trigger for dark radiation free streaming approximation (idm-idr) */
  double tau_idr;                /**< decoupling tau for idr */
  double tau_idm_dr;             /**< decoupling tau for idm_dr */
  //@}

  /** @name - initial conformal time at which thermodynamical variables have been be integrated */

  //@{

  double tau_ini; /**< initial conformal time at which thermodynamical variables have been be integrated */

  //@}

  /** @name - other thermodynamical quantities */

  //@{

  double fHe;  /**< \f$ f_{He} \f$: primordial helium-to-hydrogen nucleon ratio 4*n_He/n_H */
  double n_e;  /**< total number density of electrons today (free or not) */

  //@}

 /** START #TWIN SECTOR */
  double n_e_twin; /**< total number density of twin electrons today (free or not) */
  double z_rec_twin;   /**< z at which the visibility reaches its maximum (= recombination redshift) */
  double tau_rec_twin; /**< conformal time at which the visibility reaches its maximum (= recombination time) */
  double rs_rec_twin;  /**< comoving sound horizon at recombination */
  double ds_rec_twin;  /**< physical sound horizon at recombination */
  double ra_rec_twin;  /**< conformal angular diameter distance to recombination */
  double da_rec_twin;  /**< physical angular diameter distance to recombination */

  double tau_d_twin; /**< twin baryon drag time **/
  double z_d_twin;  /**<twin baryon drag redshift **/
  double ds_d_twin; /**< physical sound horizon at twin baryon drag */
  double rs_d_twin; /**< comoving sound horizon at twin baryond drag */
    
  double angular_rescaling_twin; /**< [ratio ra_rec / (tau0-tau_rec)]: gives CMB rescaling in angular space relative to flat model (=1 for curvature K=0) */
  double tau_free_streaming_twin;     /**< minimum value of tau at which free-streaming approximation can be switched on */
  double Z_H_REC_MAX_twin; /**< Redshift after which twin H recombination happens */
  int index_tau_twin; /* index of tau_free_streaming_twin */
  /** END TWIN SECTOR */


  /**
   *@name - some flags needed for thermodynamics functions
   */

  //@{

  short inter_normal;  /**< flag for calling thermodynamics_at_z and find position in interpolation table normally */
  short inter_closeby; /**< flag for calling thermodynamics_at_z and find position in interpolation table starting from previous position in previous call */

  //@}

  /** @name - technical parameters */

  //@{

  short thermodynamics_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */
  short hyrec_verbose; /**< flag regulating the amount of information sent to standard output from hyrec (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}

};

/**
 * Other structures that are used during the thermodynamics module
 * execution (i.e. during thermodynamics_init()) but get erased later
 * on: thus they cannot be accessed by other modules.
 */

/**
 * Vector of thermodynamical quantities to integrate over, and indices of this vector
 */



struct thermo_vector {

  int ti_size;          /**< size of thermo vector (ti stands for thermodynamical, integrated) */

  int index_ti_x_H;     /**< index for hydrogen fraction in y */
  int index_ti_x_He;    /**< index for helium fraction in y */
  int index_ti_D_Tmat;  /**< index for T_mat - T_photon [Kelvin] in y */
  /* BEGIN #TWIN SECTOR */
  int index_ti_x_H_twin;     /**< index for twin hydrogen fraction in y */
  int index_ti_x_He_twin;    /**< index for twin helium fraction in y */
  int index_ti_D_Tmat_twin;  /**< index for T_mat_twin - T_photon_twin [Kelvin] in y */
  /* END TWIN SECTOR */ 

  double * y;           /**< vector of quantities to be integrated */
  double * dy;          /**< time-derivative of the same vector */

  int * used_in_output; /**< boolean array specifying which quantities enter in the calculation of output functions */

};

/**
 * Workspace for differential equation of thermodynamics
 */

struct thermo_diffeq_workspace {

  double x_H;        /**< Hydrogen ionization fraction */
  double x_He;       /**< Helium ionization fraction */
  double x_noreio;   /**< Electron ionization fraction, not taking into account reionization */
  double x_reio;     /**< Electron ionization fraction, taking into account reionization */

  double x;          /**< total ionization fraction following usual CMB convention, n_free/n_H = x_H + fHe * x_He; */

  double Tmat;       /**< matter temperature */

  /* index of approximation schemes for the thermal history */
  int index_ap_brec; /**< before all recombination */
  int index_ap_He1;  /**< during 1st He-recombination (HeIII) */
  int index_ap_He1f; /**< in between 1st and 2nd He recombination */
  int index_ap_He2;  /**< beginning of 2nd He-recombination (HeII) */
  int index_ap_H;    /**< beginning of H-recombination (HI) */
  int index_ap_frec; /**< during and after full H- and HeII-recombination */
  int index_ap_reio; /**< during reionization */

  int ap_current;     /** current approximation scheme index */
  int ap_size;        /**< number of approximation intervals used during evolver loop */
  int ap_size_loaded; /**< number of all approximations  */



  double * ap_z_limits;       /**< vector storing ending limits of each approximation */
  double * ap_z_limits_delta; /**< vector storing smoothing deltas of each approximation */

  int require_H;  /** in given approximation scheme, do we need to integrate hydrogen ionization fraction? */
  int require_He; /** in given approximation scheme, do we need to integrate helium ionization fraction? */
   
    
  /* BEGIN #TWIN SECTOR */  
  double x_H_twin;        /**< Twin hydrogen ionization fraction */
  double x_He_twin;       /**< Twin helium ionization fraction */
  double x_noreio_twin;   /**< Twin electron ionization fraction, not taking into account reionization - no twin reionization */

  double x_twin;          /**< total twin ionization fraction following usual CMB convention, n_free_twin/n_H_twin = x_H_twin + fHe_twin * x_He_twin; */
  double x_2s_twin;       /* Fraction of twin electrons in 2s state */
  double Tmat_twin;       /**< Twin matter temperature */
  /* Additional indices of approximation schemes for the thermal history including the dark sector */
  int index_ap_brec_twin;
  int index_ap_He1_twin;  /**< during 1st twin He-recombination (HeIII) */
  //int index_ap_He1f_twin; /**< in between 1st and 2nd twin He recombination */
  int index_ap_He2_twin;  /**< beginning of 2nd twin He-recombination (HeII) */
  int index_ap_H_twin;    /**< beginning of twin H-recombination (HI) */
  int index_ap_frec_twin; /**< during and after full twin H- and twin HeII-recombination, up to beginning of SM helium recombination */
    
  int ap_current_twin;     /** current approximation scheme index */
  int ap_size_twin;        /**< number of approximation intervals used during evolver loop */
  int ap_size_loaded_twin; /**< number of all approximations  */
  double * ap_z_limits_twin;       /**< vector storing ending limits of each approximation */
  double * ap_z_limits_delta_twin; /**< vector storing smoothing deltas of each approximation */    
  int require_H_twin;  /** in given approximation scheme, do we need to integrate hydrogen ionization fraction? */
  int require_He_twin; /** in given approximation scheme, do we need to integrate helium ionization fraction? */
  double Fpr_tab[200][40]; /* Table of values of photorecombination atomic physics function */
  double Fpi_tab[200]; /* List of values of photoionization atomic physics function */
  double ToB_tab[200]; /* List of values of TR/B_D that F_pi and F_pr are evaluated at*/
  double TMoTR_tab[40]; /* List of values of T_m/T_R that F_pr is evaluated at */
  double A2s_tab[200]; /* List of values of photorecombination coefficient for 2s state, with Tm=Tr */
  double A2s_ToB_tab[200]; /* List of values of T_R/B_D that A2s is evaluated at */
  /* END TWIN SECTOR */
    
  struct thermo_vector * ptv;       /**< pointer to vector of integrated quantities and their time-derivatives */
  struct thermohyrec * phyrec;     /**< pointer to wrapper of HyRec structure */
  struct thermohyrec * phyrec_twin; /**< pointer to wrapper of HyRec structure for twin sector */
  struct thermorecfast * precfast; /**< pointer to wrapper of RecFast structure */

};

/**
 * Workspace for reionization
 */

struct thermo_reionization_parameters{

  /* parameters used by reio_camb */

  int index_re_reio_redshift;  /**< hydrogen reionization redshift */
  int index_re_reio_exponent;  /**< an exponent used in the function x_e(z) in the reio_camb scheme */
  int index_re_reio_width;     /**< a width defining the duration of hydrogen reionization in the reio_camb scheme */
  int index_re_xe_before; /**< ionization fraction at redshift 'reio_start' */
  int index_re_xe_after;  /**< ionization fraction after full reionization */
  int index_re_helium_fullreio_fraction; /**< helium full reionization fraction inferred from primordial helium fraction */
  int index_re_helium_fullreio_redshift; /**< helium full reionization redshift */
  int index_re_helium_fullreio_width;    /**< a width defining the duration of helium full reionization in the reio_camb scheme */

  /* parameters used by reio_bins_tanh, reio_many_tanh, reio_inter */

  int re_z_size;                /**< number of reionization jumps */
  int index_re_first_z;        /**< redshift at which we start to impose reionization function */
  int index_re_first_xe;       /**< ionization fraction at redshift first_z (inferred from recombination code) */
  int index_re_step_sharpness; /**< sharpness of tanh jump */

  /* parameters used by all schemes */

  int index_re_reio_start;     /**< redshift above which hydrogen reionization neglected */

  double * reionization_parameters; /**< vector containing all reionization parameters necessary to compute xe(z) */
  int re_size;              /**< length of vector reionization_parameters */
};

/**
 * General parameters relevant to thermal history and pointers to few other more specialised worspaces
 */

struct thermo_workspace {

  /* Number of z values */
  int Nz_reco_lin;             /**< number of redshifts linearly sampled for recombination during the evolver loop */
  int Nz_reco_log;             /**< number of redshifts logarithmically sampled for recombination during the evolver loop */
  int Nz_reco;                 /**< number of redshifts for recombination during the evolver loop */
  int Nz_reio;                 /**< number of redshift points of reionization during evolver loop*/
  int Nz_tot;                  /**< total number of sampled redshifts */

  /* Most important and useful parameters of evolution */
  double YHe;          /**< defined as in RECFAST : primordial helium mass fraction */
  double fHe;          /**< defined as in RECFAST : primordial helium-to-hydrogen nucleon ratio */
  double SIunit_H0;    /**< defined as in RECFAST : Hubble parameter today in SI units */
  double SIunit_nH0;   /**< defined as in RECFAST : Hydrogen number density today in SI units*/
  double Tcmb;         /**< CMB temperature today in Kelvin */

  /* Most important and useful constants */
  double const_NR_numberdens;  /**< prefactor in number density of nonrelativistic species */
  double const_Tion_H;         /**< ionization energy for HI as temperature */
  double const_Tion_HeI;       /**< ionization energy for HeI as temperature */
  double const_Tion_HeII;      /**< ionization energy for HeII as temperature */

  double reionization_optical_depth; /**< reionization optical depth inferred from reionization history */

  int last_index_back; /**< nearest location in background table */
    
  /* BEGIN #TWIN SECTOR */
  double YHe_twin;
  double fHe_twin;
  double SIunit_H0_twin;
  double SIunit_nH0_twin;
  double Tnow_twin;
  /* Most important and useful constants */
  double const_NR_numberdens_twin;  /**< prefactor in number density of nonrelativistic twin species */
  double const_Tion_H_twin;         /**< ionization energy for twin HI as temperature */
  double const_Tion_HeI_twin;       /**< ionization energy for twin HeI as temperature */
  double const_Tion_HeII_twin;      /**< ionization energy for twin HeII as temperature */
  double z_decoupling_estimate;  /** Estimate of the redshift of dark decoupling, used to define approx switch. **/
  double z_H_twin_boltzmann_trigger; /**< Redshift where twin hydrogen should start being evolved with Boltzmann equation **/
  double z_H_twin_saha_trigger;      /**< Redshift where twin hydrogen should start being evolved with Sahaequation **/
  double z_He1_twin_trigger;         /**< Redshift where twin HeII->HeI should start being evolved with Saha equation **/
  //double z_He1f_twin_trigger;
  double z_He2_twin_trigger;         /**< Redshift where twin HeIII->HeII should start being evolved with Saha equation **/
  double RLya_twin_min;
  double dxdt_at_RLya_min;
  double dxdt_max;
    
    
  double Tmat_decoupled_reference_twin;
  double z_decoupled_reference_twin;
  /* END TWIN SECTOR */
    
  /* START #TWIN SECTOR */
  short nodarkrecomb_twin;  
    
  /* END TWIN SECTOR */  
  
  struct thermo_diffeq_workspace * ptdw;        /**< pointer to workspace for differential equations */
  struct thermo_reionization_parameters * ptrp; /**< pointer to workspace for reionization */

};

/**
 * temporary parameters and workspace passed to the thermodynamics_derivs function
 */

struct thermodynamics_parameters_and_workspace {

  /* structures containing fixed input parameters (indices, ...) */
  struct background * pba;
  struct precision * ppr;
  struct thermodynamics * pth;

  /* workspace */
  struct thermo_workspace * ptw;
  double * pvecback;
};

/**************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  /* external functions of the module*/

  int thermodynamics_at_z(struct background * pba,
                          struct thermodynamics * pth,
                          double z,
                          enum interpolation_method inter_mode,
                          int * last_index,
                          double * pvecback,
                          double * pvecthermo);

  int thermodynamics_init(struct precision * ppr,
                          struct background * pba,
                          struct thermodynamics * pth);

  int thermodynamics_free(struct thermodynamics * pth);

  /* internal functions of the module */

  int thermodynamics_helium_from_bbn(struct precision * ppr,
                                     struct background * pba,
                                     struct thermodynamics * pth);

  int thermodynamics_checks(struct precision * ppr,
                            struct background* pba,
                            struct thermodynamics * pth);

  int thermodynamics_workspace_init(struct precision * ppr,
                                    struct background * pba,
                                    struct thermodynamics * pth,
                                    struct thermo_workspace * ptw);

  int thermodynamics_indices(struct background * pba,
                             struct thermodynamics * pth,
                             struct thermo_workspace* ptw);

  int thermodynamics_lists(struct precision * ppr,
                           struct background* pba,
                           struct thermodynamics* pth,
                           struct thermo_workspace* ptw);

  int thermodynamics_set_parameters_reionization(struct precision * ppr,
                                                 struct background * pba,
                                                 struct thermodynamics * pth,
                                                 struct thermo_reionization_parameters * preio);

  int thermodynamics_solve(struct precision * ppr,
                           struct background * pba,
                           struct thermodynamics * pth,
                           struct thermo_workspace* ptw,
                           double * pvecback);

  int thermodynamics_calculate_remaining_quantities(struct precision * ppr,
                                                    struct background * pba,
                                                    struct thermodynamics* pth,
                                                    double* pvecback);

  int thermodynamics_output_summary(struct background* pba,
                                    struct thermodynamics* pth);

  int thermodynamics_workspace_free(struct thermodynamics* pth, struct thermo_workspace * ptw);

  int thermodynamics_vector_init(struct precision * ppr,
                                 struct background * pba,
                                 struct thermodynamics * pth,
                                 double z,
                                 struct thermo_workspace * ptw);

  int thermodynamics_reionization_evolve_with_tau(struct thermodynamics_parameters_and_workspace * tpaw,
                                                  double mz_ini,
                                                  double mz_end,
                                                  double * mz_output,
                                                  int Nz);

  int thermodynamics_derivs(
                            double mz,
                            double * y,
                            double * dy,
                            void * parameters_and_workspace,
                            ErrorMsg error_message
                            );

  int thermodynamics_timescale(double z,
                               void * thermo_parameters_and_workspace,
                               double * timescale,
                               ErrorMsg error_message);

  int thermodynamics_sources(double mz,
                             double * y,
                             double * dy,
                             int index_z,
                             void * thermo_parameters_and_workspace,
                             ErrorMsg error_message);

  int thermodynamics_reionization_get_tau(struct precision * ppr,
                                          struct background * pba,
                                          struct thermodynamics * pth,
                                          struct thermo_workspace * ptw);

  int thermodynamics_vector_free(struct thermo_vector * tv);


  int thermodynamics_calculate_conformal_drag_time(struct background* pba,
                                                   struct thermodynamics* pth,
                                                   double* pvecback);

  int thermodynamics_calculate_damping_scale(struct background* pba,
                                             struct thermodynamics* pth,
                                             double* pvecback);

  int thermodynamics_calculate_opticals(struct background* pba,
                                        struct precision* ppr,
                                        struct thermodynamics* pth);

  int thermodynamics_calculate_idm_dr_quantities(struct precision * ppr,
                                                 struct background * pba,
                                                 struct thermodynamics * pth,
                                                 double* pvecback);

  int thermodynamics_calculate_recombination_quantities(struct precision* ppr,
                                                        struct background * pba,
                                                        struct thermodynamics* pth,
                                                        double* pvecback);

  int thermodynamics_calculate_drag_quantities(struct precision* ppr,
                                               struct background * pba,
                                               struct thermodynamics* pth,
                                               double* pvecback);

  int thermodynamics_ionization_fractions(
                                          double z,
                                          double * y,
                                          struct background * pba,
                                          struct thermodynamics * pth,
                                          struct thermo_workspace * ptw,
                                          int current_ap
                                          );

  int thermodynamics_reionization_function(double z,
                                           struct thermodynamics * pth,
                                           struct thermo_reionization_parameters * preio,
                                           double * x);

  int thermodynamics_output_titles(struct background * pba,
                                   struct thermodynamics *pth,
                                   char titles[_MAXTITLESTRINGLENGTH_]);

  int thermodynamics_output_data(struct background * pba,
                                 struct thermodynamics *pth,
                                 int number_of_titles,
                                 double *data);


  /** START #TWIN SECTOR */
  /* Twin thermodynamics functions */
  int thermodynamics_twin_helium_from_bbn(
                   struct precision * ppr,
                   struct background * pba,
                   struct thermodynamics * pth
                   );
  int thermodynamics_sources_twin(double mz,
                             double * y,
                             double * dy,
                             int index_z,
                             void * thermo_parameters_and_workspace,
                             ErrorMsg error_message);
    
  int thermodynamics_vector_init_twin(struct precision * ppr,
                                 struct background * pba,
                                 struct thermodynamics * pth,
                                 double z,
                                 struct thermo_workspace * ptw);
    
  int thermodynamics_derivs_twin(
                              double mz,
                              double * y,
                              double * dy,
                              void * parameters_and_workspace,
                              ErrorMsg error_message
                              );
  int thermodynamics_ionization_fractions_twin(
                                          double z,
                                          double * y,
                                          struct background * pba,
                                          struct thermodynamics * pth,
                                          struct thermo_workspace * ptw,
                                          int current_ap_twin
                                          );
  double my_1D_interp(double x, double * xi, double * yi);
  double my_2D_interp(double x, double y, double * xi, double * yi, double zi[200][40]);
  /** END TWIN SECTOR */


#ifdef __cplusplus
}
#endif

/**************************************************************/

/**
 * @name some flags
 */

//@{

#define _YHE_BBN_ -1 /**< value assigned to the parameter pth->YHe by the input module when this parameter must be computed using BBN tables */

//@}

/**
 * @name Some basic constants needed by RECFAST:
 */

//@{

#define _m_e_ 9.10938215e-31  /**< electron mass in Kg */
#define _m_p_ 1.672621637e-27 /**< proton mass in Kg */
#define _m_H_ 1.673575e-27    /**< Hydrogen mass in Kg */
#define _not4_ 3.9715         /**< Helium to Hydrogen mass ratio */
#define _sigma_ 6.6524616e-29 /**< Thomson cross-section in m^2 */

#define _RECFAST_INTEG_SIZE_ 3

//@}

/**
 * @name Some other physical constants
 */

//@{

#define _s_over_Mpc_ 9.71561189e-15  /**< conversion factor from s to megaparsecs (1 s= const*Mpc) */
#define _Mpc_over_GeV_ 1.56373832e38  /**< conversion factor from GeV to megaparsecs (1 GeV= const/Mpc) */
#define _GeV_over_kg_ 1.78266191e-27  /**< conversion factor from GeV to kg  (1 GeV= const*kg) */
#define _GeVcm3_over_Mpc2_ 94.7024726  /**< conversion factor from  CLASS_rho 1/Mpc^2 to rho in GeV/cm^3 (rho in GeV/cm^3=const*CLASS_rho) */
#define _Jm3_over_Mpc2_ 0.0151730087  /**< conversion factor from  CLASS_rho 1/Mpc^2 to rho in Joule/m^3 (rho in Joule/m^3=const*CLASS_rho) */
#define _Sun_mass_ 1.98855e30 /**< sun mass in kg */
#define _eV_over_Kelvin_ 8.61733034e-5   /**< kB in eV/K */
#define _eV_over_joules_ 6.24150647996e+18 /**< eV/J */


//@}

/* @endcond */

/**
 * @name Some limits imposed on cosmological parameter values:
 */

//@{

#define _YHE_BIG_ 0.5      /**< maximal \f$ Y_{He} \f$ */
#define _YHE_SMALL_ 0.01   /**< minimal \f$ Y_{He} \f$ */
#define _Z_REC_MAX_ 2000.
#define _Z_REC_MIN_ 500.

//@}

/** START #TWIN SECTOR: April 26 redefined sigma-twin and m_H_twin for adm. */
#define _mu_twin_ (pba->m_e_dark*_GeV_over_kg_*pba->m_p_dark*_GeV_over_kg_)/((pba->m_e_dark*_GeV_over_kg_) + (pba->m_p_dark*_GeV_over_kg_)) /** Dark reduced mass */
#define _m_H_twin (pba->m_e_dark*_GeV_over_kg_ + pba->m_p_dark*_GeV_over_kg_ -0.5 * _mu_twin_ * pba->alphafs_dark * pba->alphafs_dark) /** Dark H mass*/
#define _m_e_twin (_m_e_*pba->ratio_vev_twin) /** Twin electron mass*/
#define _epsilon0_perm_ 8.8541878128e-12 /** Vacuum Permittivity*/
#define _sigma_twin  (_sigma_*(pba->alphafs_dark/0.00729735)*(pba->alphafs_dark/0.00729735)/pow(pba->ratio_vev_twin,2)) /**< Twin Thomson cross-section in m^2 */
#define _Z_REC_MIN_twin 0.1 /* CHECK */
#define _L_H_ion_twin (1.096787737e7*(_mu_twin_/_m_e_)*(pba->alphafs_dark/0.00729735)*(pba->alphafs_dark/0.00729735))
#define _L_He1_ion_twin (1.98310772e7*pba->ratio_vev_twin*(pba->alphafs_dark/0.00729735)*(pba->alphafs_dark/0.00729735))
#define _L_He2_ion_twin (4.389088863e7*pba->ratio_vev_twin*(pba->alphafs_dark/0.00729735)*(pba->alphafs_dark/0.00729735))
/** END TWIN SECTOR */

#endif
