/////////////////////////////////////////////////////
// Antarctic krill assessment model
// NOAA/SWFSC/AERD
// December 4, 2017
//
//  merging catch_opt options from Aug 9 tpl
//  with converging tpl code from Jul 29
//  1 = Ianelli F; 2 = subtract catch; 3 = Pope's approximation;
//  4 = Thomson/Bravington F
//
//  added to Fmort_pen:
//    if (active(log_avg_fmort))   
//  so that fpen won't be calculated when fmort isn't being calculated
//
//  put non-DATA_SECTION couts in if DebugOut > 0 statement to speed up mcmc
//  srv_like is pre-Francis (2011) January equation
//  M penalty ON in compute_priors
//  all estimated parameters in MCMC ouput 
//  added 2 forms of Ricker recruitment, fixed rec_pen cout error
//  Compute_priors for sigmar, steepness
//  added Bzero, Btot to sdreport_matrix


  // ============================
DATA_SECTION
  // ============================

  !!CLASS ofstream mceval("mceval.dat")

  int i;
  int j;
  int iyr;
  int isrv;
  int ifsh;
  int nfsh;
  int age;
  int ilen;
  int iarea;
  int jarea;
  int istart;
  int imov;
  init_int      DebugOut;
  init_int      Terminal_phase;                   // Last phase for the model
  init_int      styr;                             // First model year
  init_int      endyr;                            // Last model year
  init_int      nproj_yrs;
  init_int      mat_yrs                           // Number of years with maturity data
  init_int      rec_age;
  init_int      oldest_age;                       // Last age-class
  init_int      nl_bins;                          // Number of length bins -1
  init_ivector  lbins(1,nl_bins);                 // Length-bins for compositions
  int           end_datyrs;
  int           nal_bins;
  !!            nal_bins = nl_bins+1;
  int           nages;                            // Age-range
  !!            nages = oldest_age-rec_age+1;     // Number of ages
  int           iters;                            // model iterations
  int           styr_rec;
  int           styr_sp;
  int           endyr_sp;
  int           nyrs;
  int           ndatyrs;
  int           nmov;
  !! nyrs = endyr-styr+1;  
  !! end_datyrs = endyr-nproj_yrs;
  !! ndatyrs = end_datyrs-styr+1;
  !! iters = 0;
  !! styr_rec = (styr - nages) + 1;               // First year of recruitment
  !! styr_sp  = styr_rec-rec_age;                 // First year of spawning biomass  
  !! endyr_sp = endyr-rec_age  - 1;               // endyr year of (main) spawning biomass
  init_int      nareas;                           // Number of areas
  !! if (nareas > 1)
  !!   nmov = nareas*(nareas-1);                  // number of movement rates to estimate
  !! if (nareas == 1)
  !!   nmov = nareas;                             // number of movement rates to estimate
  init_int      nsrv;                             // Number of surveys
  init_int      nsel;                             // Number of different selectivity patterns
  init_vector   msrv_id(1,nsrv);                   // Identifier for the selectivity pattern for each survey
  int           msrv;                             // index to hold selectivity pattern in nsrv loop
  init_vector   mo_srv_leg(1,nsrv);               // Months in which survey-legs occur
  vector        rho(1,nsrv);
  !! for(isrv=1;isrv<=nsrv;isrv++)
  !!        rho(isrv) = (mo_srv_leg(isrv)-1)/12;
  !! cout << "rho = " << rho << endl;
  !! cout << "styr = " << styr << endl;
  !! cout << "endyr = " << endyr << endl;
  !! cout << "styr_rec = " << styr_rec << endl;
  !! cout << "styr_sp = " << styr_sp << endl;
  !! cout << "nages = " << nages << endl;
  !! cout << "nsrv = " << nsrv << endl;
  !! cout << "nsel = " << nsel << endl;
  !! cout << "msrv_id = " << msrv_id << endl;
  !! cout << "nareas = " << nareas << endl;
  !! cout << "mo_srv_leg = " << mo_srv_leg << endl;
  !! cout << "nmov = " << nmov << endl;
  number Steepness_UB;                            // Upper bound for Steepness
  !! Steepness_UB =   .999;                               
  vector age_vector(1,nages);
    !! for (age=1;age<=nages;age++)
      !!  age_vector(age) = double(age);

  // ============================
  //Read (net tow) survey data on length compositions
  init_imatrix nyrs_srv_comp(1,nsrv,1,nareas)            // Number of years of composition data for surveys
  !! cout <<   "nyrs_srv_comp = " << endl << nyrs_srv_comp << endl;
  !! int lnsrv = nsrv;
  !! int lnareas = nareas;
  !! imatrix lnyrs_srv_comp = nyrs_srv_comp;
  init_3darray yrs_srv_comp(1,lnsrv,1,lnareas,1,lnyrs_srv_comp);    // Years with survey compostions
  !! cout <<   "yrs_srv_comp = " << endl << yrs_srv_comp << endl;
  init_3darray nsmpl_srv(1,lnsrv,1,lnareas,1,lnyrs_srv_comp);// Effective sample sizes of composition survey years
  !! cout <<   "nsmpl_srv = " << endl << nsmpl_srv << endl;
  init_vector  wt_len(1,nl_bins)
  init_vector  wt_age_init(1,nages)
  init_vector  mat_age(1,nages)
  !! int             lnl_bins=nl_bins;
  init_4darray oc_srv(1,lnsrv,1,lnareas,1,lnyrs_srv_comp,1,lnl_bins);// Survey length compositions
  !! cout <<   "wt_len = " << endl << wt_len << endl;
  !! cout <<   "wt_age_init = " << endl << wt_age_init << endl;
  !! cout <<   "mat_age = " << endl << mat_age << endl;
  !! for(isrv=1;isrv<=nsrv;isrv++)
  !! for(iarea=1;iarea<=nareas;iarea++)
  !!   cout << "oc_srv(" << isrv << ", "<< iarea << ") = " << endl << oc_srv(isrv,iarea) << endl;

  // ============================
  //Read maturities at length
  init_3darray mat_len(1,1,1,nl_bins,1,mat_yrs);
  !! cout <<   "mat_len = " << endl << mat_len << endl;
  vector       avg_mature(1,nl_bins);
  !!for (iarea=1; iarea<=1; iarea++)
    !!{
      !!for (ilen=1; ilen<=nl_bins; ilen++)
      !!{
      !! avg_mature(ilen) = mean(mat_len(iarea,ilen));
      !!}
    !!}
  !! cout << endl << "avg_mature = " << endl << avg_mature << endl;

  // =========================                  
  // Read survey-index data
  init_imatrix  nyrs_srv(1,nsrv,1,nareas);                      // Number of years of biomass for each survey
  !! imatrix lnyrs_srv = nyrs_srv;                              // ragged array magic
  init_3darray  yrs_srv(1,lnsrv,1,lnareas,1,lnyrs_srv)          // Years with biomass survey data
  init_3darray   obs_srv(1,lnsrv,1,lnareas,1,lnyrs_srv)         // The biomass survey indices
  init_3darray   obs_se_srv(1,lnsrv,1,lnareas,1,lnyrs_srv)      // The biomass survey SDs
  !! cout <<   "nyrs_srv = " << endl << nyrs_srv << endl;
  !! cout <<   "yrs_srv = " << endl << yrs_srv << endl;
  !! cout <<   "obs_srv = " << endl << obs_srv << endl;
  !! cout <<   "obs_se_srv = " << endl << obs_se_srv << endl;

  // =========================                  
  // Read fisheries data
  // =========================                  
  init_vector    nfsh_area(1,nareas);                            // Number of fisheries (fleets) per area
  !! nfsh = sum(nfsh_area);                                     // Total number of fleets, all area combinations
  init_vector    areas_fsh(1,nfsh);                             // area for each fleet
  init_vector    mo_fsh_leg(1,nfsh);                            // Months in which survey-legs occur
  vector         tau(1,nfsh);
  !! for(ifsh=1;ifsh<=nfsh;ifsh++)
  !!             tau(ifsh) = (mo_fsh_leg(ifsh)-1)/12;
  init_ivector   nyrs_fsh_comp(1,nfsh);                         // Number of years of composition data for fisheries
  !! cout <<    "nyrs_fsh_comp = " << endl << nyrs_fsh_comp << endl;
  //!! int lnfsh = nfsh;
  //!! imatrix     lnyrs_fsh_comp = nyrs_fsh_comp;
  init_matrix    yrs_fsh_comp(1,nfsh,1,nyrs_fsh_comp);        // Years with survey compostions
  init_matrix    nsmpl_fsh(1,nfsh,1,nyrs_fsh_comp);           // Effective sample sizes of fishery composition years
  init_3darray   oc_fsh(1,nfsh,1,nyrs_fsh_comp,1,nl_bins+13);// fisheries length compositions -dhk: fix l_bin index
  !! for(ifsh=1;ifsh<=nfsh;ifsh++)
  !!   cout << "oc_fsh(" << ifsh << ") = " << endl << oc_fsh(ifsh) << endl;
  // =========================                  
  // Read fisheries-index data
  init_ivector  nyrs_fsh(1,nfsh);                          // Number of years of biomass for each fisheries
  //!! imatrix lnyrs_fsh = nyrs_fsh;                                  // ragged array magic
  init_matrix    yrs_fsh(1,nfsh,1,nyrs_fsh);               // Years with biomass fisheries data
  init_matrix    obs_catch(1,nfsh,styr,endyr);             // The biomass fisheries indices
  init_number    gamma;                                    // no longer used
  init_vector    obs_se_catch(1,nareas)          // The biomass fisheries SDs
  vector         obs_catch_pen(1,nareas)                    // Convert cv_catchbiomass to penalty
  !! for(iarea=1;iarea<=nareas;iarea++)
  !!  obs_catch_pen(iarea) = 1.0/(2*obs_se_catch(iarea)*obs_se_catch(iarea));
  matrix         effort(1,nfsh,styr,endyr)
  !! for(ifsh=1;ifsh<=nfsh;ifsh++)
  !!   for(iyr=styr;iyr<=endyr;iyr++)
  !!     if(obs_catch(ifsh,iyr) > 0)
  !!       effort(ifsh,iyr) = 1;
  !!     else effort(ifsh,iyr) = 0;
  !! cout <<   "nyrs_fsh = " << endl << nyrs_fsh << endl;
  !! cout <<   "yrs_fsh = " << endl << yrs_fsh << endl;
  !! cout <<   "obs_catch = " << endl << obs_catch << endl;
  !! cout <<   "obs_se_catch = " << endl << obs_se_catch << endl;
  !! cout <<   "effort = " << endl << effort << endl;

  // ============================
  // Modeling options (.CTR values)
  // ============================

  init_int    styr_rec_est
  init_int    endyr_rec_est
  int         nrecs_est;
  !! nrecs_est = endyr_rec_est-styr_rec_est+1;
  !! cout<<"Rec estimated in styr_rec endyr: " << styr_rec    <<" "<< endyr        <<" "<<endl;
  !! cout<<" SR Curve fit  in styr_rec_est endyr: " << styr_rec_est <<" "<< endyr_rec_est <<" "<<endl;
  !! cout<<"            Model styr endyr: " <<styr        <<" "<<endyr        <<" "<<endl;
  !! cout << endl;
  init_number  natmortprior;
  init_number  cvnatmortprior;
  !! cout << "natmortprior = " << endl << natmortprior << endl;
  !! cout << "cvnatmortprior = " << cvnatmortprior << endl;


  init_vector  log_mov_init(1,nmov);
  init_matrix  rec_dev_init(1,nareas,styr_rec,endyr) 
  init_number  mean_log_rec_init;
  init_number  log_Rzero_init;
  init_number  log_influx_init;
  init_number  Mdevs_init;

  init_int    use_priors;                             // 1 if Compute_priors
  init_vector qprior(1,nsrv);                         // Prior for survey-q
  vector      log_qprior(1,nsrv);
  !! log_qprior = log(qprior);
  init_vector cvqprior(1,nsrv);
  !! cout << "rec_dev_init = " << endl << rec_dev_init << endl;
  !! cout << "log_mov_init = " << log_mov_init << endl;
  !! cout << "mean_log_rec_init = " << mean_log_rec_init << endl;
  !! cout << "log_Rzero_init = " << log_Rzero_init << endl;
  !! cout << "log_influx_init = " << log_influx_init << endl;
  !! cout << "Mdevs_init = " << Mdevs_init << endl;
  !! cout << "qprior = " << qprior << endl;
  !! cout << "log_qprior = " << log_qprior << endl;
  !! cout << "cvqprior = " << cvqprior << endl;
  init_int     influx_on;                    // 1 = estimate influx, 0 = estimation off
  init_int     phase_M;
  init_int     phase_Mdevs;
  init_int     phase_Rzero;                  // Phase for Rzero
  init_int     phase_srec;
  init_int     phase_meanRec;                // Phase for mean_log_rec
  init_int     phase_recdevs;                // Phase for rec_dev
  init_int     phase_al;                     // Phase for vB age-length
  init_int     phase_move_m;                 // Phase for mean movement
  init_int     phase_move_d;                 // Phase for movement deviations
  init_int     phase_sigmar;                 // Phase for recruitment variability
  init_int     phase_fmort;                  // Phase for fishing mortality
  init_int     catch_opt;                  // Phase for catch estimation
  init_ivector phase_q_srv(1,nsrv);          // Phase for survey q
  init_ivector phase_q_fsh(1,nfsh);          // Phase for fishery q
  init_int     SrType;                       // Stock-Recruit type: 2 BHolt, 1 Ricker (Dorn), 4 Ricker(old)

  init_number steepnessprior;                // Prior for steepness
  init_number cvsteepnessprior;

  init_ivector srv_sel_opt(1,nsrv);
  init_int     q_age_min                     // Minimum age for variable survey selectivity
  init_int     q_age_max                     // Maximum age for variable survey selectivity  
  !! cout << "influx_on = " << influx_on << endl;
  !! cout << "phase_M = " << phase_M << endl;
  !! cout << "phase_move_m = " << phase_move_m << endl;
  !! cout << "phase_Rzero = " << phase_Rzero << endl;
  !! cout << "phase_srec = " << phase_srec << endl;
  !! cout << "phase_meanRec = " << phase_meanRec << endl;
  !! cout << "phase_recdevs = " << phase_recdevs << endl;
  !! cout << "phase_al = " << phase_al << endl;
  !! cout << "phase_move_m = " << phase_move_m << endl;
  !! cout << "phase_fmort = " << phase_fmort << endl;
  !! cout << "phase_sigmar = " << phase_sigmar << endl;
  !! cout << "phase_q_srv = " << phase_q_srv << endl;
  !! cout << "phase_q_fsh = " << phase_q_fsh << endl;
  !! cout <<  "SrType = " << SrType << endl;
  !! cout <<   "steepnessprior = " << steepnessprior << endl;
  !! cout <<   "cvsteepnessprior = " << cvsteepnessprior << endl;
  !! cout <<   "steepnessprior = " << steepnessprior << endl;
  !! cout <<   "srv_sel_opt = " << srv_sel_opt << endl;
  !! cout <<   "q_age_min, q_age_max =" << endl;
  !! cout <<   q_age_min << " " << q_age_max << endl;
  ivector      nselages_in_srv(1,nsrv);      //  Number of age classes with selectivities by survey
  vector       offset_srv(1,nsrv); 
  int          seldecage ;
  !! seldecage = int(nages/2);
  int          phase_sel_srv
  vector       curv_pen_srv(1,nsrv)
  vector       sel_slp_in_srv(1,nsrv)
  vector       logsel_slp_in_srv(1,nsrv)
  vector       sel_inf_in_srv(1,nsrv)
  vector       logsel_inf_in_srv(1,nsrv)  // used in srv selectivity case 4 only
  vector       sel_dslp_in_srv(1,nsrv)
  vector       logsel_dslp_in_srv(1,nsrv)
  vector       sel_dinf_in_srv(1,nsrv)
  vector       seldec_pen_srv(1,nsrv) ;
  matrix       sel_change_in_srv(1,nsrv,styr,endyr);

  // Phase of estimation
  ivector      phase_selcoff_srv(1,nsrv);
  ivector      phase_logist_srv(1,nsrv);
  ivector      phase_dlogist_srv(1,nsrv);
  init_vector sigmarprior(1,nareas);          // Prior for sigma-R
  vector log_sigmarprior(1,nareas);
  !! log_sigmarprior = log(sigmarprior);
  init_vector cvsigmarprior(1,nareas);
  !! cout <<   "log_sigmarprior = " << log_sigmarprior << endl;

 LOCAL_CALCS
  sel_change_in_srv.initialize()   ;  
  for (isrv=1; isrv<=nsrv;isrv++)
  {
    msrv = msrv_id(isrv);                     // merged survey identifier
    nselages_in_srv(msrv) = nages-1;
    cout <<   "nselages_in_srv(msrv) = " << nselages_in_srv(msrv) << endl;
    cout <<   "srv_sel_opt(isrv) = " << srv_sel_opt(isrv) << endl;
    switch (srv_sel_opt(isrv))
    {
      case 1 : // Selectivity coefficients 
      {
        *(ad_comm::global_datafile) >> nselages_in_srv(msrv)   ;  
        *(ad_comm::global_datafile) >> phase_sel_srv;  
        *(ad_comm::global_datafile) >> curv_pen_srv(msrv) ;
        *(ad_comm::global_datafile) >> seldec_pen_srv(msrv) ;
        seldec_pen_srv(msrv) *= seldec_pen_srv(msrv) ;
        cout << "nselages_in_srv in area " << msrv << " = " << nselages_in_srv(msrv) << endl;
        cout << "phase_sel_srv in survey " << msrv << " = " << phase_sel_srv << endl;
        cout << "initial curv_pen_srv in area " << msrv << " = " << curv_pen_srv(msrv) << endl;
        cout << "seldec_pen_srv in area " << msrv << " = " << seldec_pen_srv(msrv) << endl;
        for (int iyr=styr;iyr<=endyr;iyr++)
          *(ad_comm::global_datafile) >> sel_change_in_srv(msrv,iyr) ;
        phase_selcoff_srv(msrv) = phase_sel_srv;
        phase_logist_srv(msrv)  = -1;
        logsel_slp_in_srv = 0. ;
        sel_inf_in_srv    = 0. ;
        phase_dlogist_srv(msrv)  = -1;
        logsel_dslp_in_srv(msrv) = 0. ;
        sel_dinf_in_srv(msrv)    = 0. ;
      }
      break;
      case 2 : // Single logistic
      {
        *(ad_comm::global_datafile) >> phase_sel_srv;  
        *(ad_comm::global_datafile) >> sel_slp_in_srv(msrv) ;
        *(ad_comm::global_datafile) >> sel_inf_in_srv(msrv) ;
        for (int iyr=styr;iyr<=endyr;iyr++) { 
          *(ad_comm::global_datafile) >> sel_change_in_srv(msrv,iyr) ; }
        phase_selcoff_srv(msrv) = -1;
        phase_logist_srv(msrv) = phase_sel_srv;
        logsel_slp_in_srv(isrv) = log(sel_slp_in_srv(isrv) +1.0e-8) ;
        phase_dlogist_srv(msrv)  = -1;
        logsel_dslp_in_srv(msrv) = 0. ;
        sel_dinf_in_srv(msrv)    = 0. ;        
        cout << "sel_slp_in_srv(" << msrv << ") = " << sel_slp_in_srv(msrv) << endl;
        cout << "sel_inf_in_srv(" << msrv << ") = " << sel_inf_in_srv(msrv) << endl;
      }
      break;
      case 3 : // Double logistic 
      {
        *(ad_comm::global_datafile) >> phase_sel_srv;  
        *(ad_comm::global_datafile) >> sel_slp_in_srv(msrv) ;
        *(ad_comm::global_datafile) >> sel_inf_in_srv(msrv) ;
        *(ad_comm::global_datafile) >> sel_dslp_in_srv(msrv) ;
        *(ad_comm::global_datafile) >> sel_dinf_in_srv(msrv) ;
        for (int iyr=styr;iyr<=endyr;iyr++) { 
          *(ad_comm::global_datafile) >> sel_change_in_srv(msrv,iyr) ; }

        phase_selcoff_srv(msrv) = -1;
        phase_logist_srv(msrv)  = phase_sel_srv;
        phase_dlogist_srv(msrv) = phase_sel_srv+1;

        logsel_slp_in_srv(isrv) = log(sel_slp_in_srv(isrv)) ;
        logsel_dslp_in_srv(msrv) = log(sel_dslp_in_srv(msrv)) ;
        cout << "sel_slp_in_srv(" << msrv << ") = " << sel_slp_in_srv(msrv) << endl;
        cout << "sel_inf_in_srv(" << msrv << ") = " << sel_inf_in_srv(msrv) << endl;
        cout << "sel_dslp_in_srv(" << msrv << ") = " << sel_slp_in_srv(msrv) << endl;
        cout << "sel_dinf_in_srv(" << msrv << ") = " << sel_inf_in_srv(msrv) << endl;
        break;
      }
      case 4 : // Double normal 
      {
        *(ad_comm::global_datafile) >> phase_sel_srv;  
        *(ad_comm::global_datafile) >> sel_slp_in_srv(msrv) ;
        *(ad_comm::global_datafile) >> sel_inf_in_srv(msrv) ;
        *(ad_comm::global_datafile) >> sel_dslp_in_srv(msrv) ;
        //*(ad_comm::global_datafile) >> sel_dinf_in_srv(msrv) ;
        for (int iyr=styr;iyr<=endyr;iyr++) { 
          *(ad_comm::global_datafile) >> sel_change_in_srv(msrv,iyr) ; }

        phase_selcoff_srv(msrv) = -1;
        phase_logist_srv(msrv)  = phase_sel_srv;
        phase_dlogist_srv(msrv) = phase_sel_srv+1;

        logsel_slp_in_srv(isrv) = log(sel_slp_in_srv(isrv)) ;
        logsel_dslp_in_srv(msrv) = log(sel_dslp_in_srv(msrv)) ;
        //logsel_inf_in_srv(msrv) = log(sel_inf_in_srv(msrv)) ; // sel50_srv for case 4
        cout << "sel_slp_in_srv(" << msrv << ") = " << sel_slp_in_srv(msrv) << endl;
        cout << "sel_inf_in_srv(" << msrv << ") = " << sel_inf_in_srv(msrv) << endl;
        cout << "sel_dslp_in_srv(" << msrv << ") = " << sel_slp_in_srv(msrv) << endl;
        cout << "sel_dinf_in_srv(" << msrv << ") (not used) = " << sel_inf_in_srv(msrv) << endl;
        break;
      }
    }
  }
 END_CALCS
  ivector n_sel_ch_srv(1,nsrv);
  imatrix yrs_sel_ch_tsrv(1,nsrv,1,endyr-styr+1);
  //!! for (isrv=1;isrv<=nsrv;isrv++) 
  //!!   if (phase_selcoff_srv(isrv)>0) 
  //!!     curv_pen_srv(isrv) = 1./ (square(curv_pen_srv(isrv))*2);
  //!! cout << "eventual curv_pen_srv(isrv) = " << curv_pen_srv(isrv) << endl;
 LOCAL_CALCS
  for (isrv=1;isrv<=nsrv;isrv++)
  {
    msrv = msrv_id(isrv);                     // merged survey identifier
    cout << "isrv = " << isrv << endl;
    cout << "msrv = " << msrv << endl;
    sel_change_in_srv(isrv,styr)=1.;
    n_sel_ch_srv(isrv)=0.;
    int j=1;
    yrs_sel_ch_tsrv(isrv,j) = styr;
    cout << "yrs_sel_ch_tsrv(" << msrv << ",j) = " << yrs_sel_ch_tsrv(msrv,j) << endl;
    for (int iyr=styr+1;iyr<=endyr;iyr++)
    {
      if(sel_change_in_srv(msrv,iyr)>0)
      {
        j++;
        yrs_sel_ch_tsrv(isrv,j) = iyr;
       cout <<"Number of select changes: survey "<<isrv<<": "<<n_sel_ch_srv(isrv)<<endl;
       cout << "yrs_sel_ch_tsrv "  << yrs_sel_ch_tsrv << endl;
      }
    }
    n_sel_ch_srv(isrv) = j;
    cout <<"Number of selectivity changes in net survey "<<isrv<<": "<<n_sel_ch_srv(isrv)<<endl;
  }
 END_CALCS
  !!cout << "sel_change_in_srv: "<< endl << sel_change_in_srv << endl;
  imatrix   yrs_sel_ch_srv(1,nsrv,1,n_sel_ch_srv);
  imatrix   nselages_srv(1,nsrv,1,n_sel_ch_srv);
 LOCAL_CALCS
  for (isrv=1; isrv<=nsrv;isrv++)
    yrs_sel_ch_srv(isrv) = yrs_sel_ch_tsrv(isrv)(1,n_sel_ch_srv(isrv));
  cout<< "yrs_sel_ch_srv: "<< endl << yrs_sel_ch_srv << endl;
 END_CALCS

 // FISHERIES SELECTIVITIES
  init_ivector fsh_sel_opt(1,nfsh);
  ivector      nselages_in_fsh(1,nfsh);      //  Number of age classes with selectivities by fishery
  vector       offset_fsh(1,nfsh); 
  int          phase_sel_fsh
  vector       curv_pen_fsh(1,nfsh)
  vector       sel_slp_in_fsh(1,nfsh)
  vector       logsel_slp_in_fsh(1,nfsh)
  vector       sel_inf_in_fsh(1,nfsh)
  vector       sel_dslp_in_fsh(1,nfsh)
  vector       logsel_dslp_in_fsh(1,nfsh)
  vector       sel_dinf_in_fsh(1,nfsh)
  vector       seldec_pen_fsh(1,nfsh) ;
  matrix       sel_change_in_fsh(1,nfsh,styr,endyr);

  // Phase of estimation
  ivector      phase_selcoff_fsh(1,nfsh);
  ivector      phase_logist_fsh(1,nfsh);
  ivector      phase_dlogist_fsh(1,nfsh);
 LOCAL_CALCS
  sel_change_in_fsh.initialize()   ;  
  for (ifsh=1; ifsh<=nfsh;ifsh++)
  {
    nselages_in_fsh(ifsh) = nages-1;
    cout <<   "nselages_in_fsh(ifsh) = " << nselages_in_fsh(ifsh) << endl;
    cout <<   "fsh_sel_opt(ifsh) = " << fsh_sel_opt(ifsh) << endl;
    switch (fsh_sel_opt(ifsh))
    {
      case 1 : // Selectivity coefficients 
      {
        *(ad_comm::global_datafile) >> nselages_in_fsh(ifsh)   ;  
        *(ad_comm::global_datafile) >> phase_sel_fsh;  
        *(ad_comm::global_datafile) >> curv_pen_fsh(ifsh) ;
        *(ad_comm::global_datafile) >> seldec_pen_fsh(ifsh) ;
        seldec_pen_fsh(ifsh) *= seldec_pen_fsh(ifsh) ;
        cout << "nselages_in_fsh in area " << ifsh << " = " << nselages_in_fsh(ifsh) << endl;
        cout << "phase_sel_fsh in survey " << ifsh << " = " << phase_sel_fsh << endl;
        cout << "initial curv_pen_fsh in area " << ifsh << " = " << curv_pen_fsh(ifsh) << endl;
        cout << "seldec_pen_fsh in area " << ifsh << " = " << seldec_pen_fsh(ifsh) << endl;
        for (int iyr=styr;iyr<=endyr;iyr++)
          *(ad_comm::global_datafile) >> sel_change_in_fsh(ifsh,iyr) ;
        phase_selcoff_fsh(ifsh) = phase_sel_fsh;
        phase_logist_fsh(ifsh)  = -1;
        logsel_slp_in_fsh = 0. ;
        sel_inf_in_fsh    = 0. ;
        phase_dlogist_fsh(ifsh)  = -1;
        logsel_dslp_in_fsh(ifsh) = 0. ;
        sel_dinf_in_fsh(ifsh)    = 0. ;
      }
      break;
      case 2 : // Single logistic
      {
        *(ad_comm::global_datafile) >> phase_sel_fsh;  
        *(ad_comm::global_datafile) >> sel_slp_in_fsh(ifsh) ;
        *(ad_comm::global_datafile) >> sel_inf_in_fsh(ifsh) ;
        for (int iyr=styr;iyr<=endyr;iyr++) { 
          *(ad_comm::global_datafile) >> sel_change_in_fsh(ifsh,iyr) ; }
        phase_selcoff_fsh(ifsh) = -1;
        phase_logist_fsh(ifsh) = phase_sel_fsh;
        logsel_slp_in_fsh(ifsh) = log(sel_slp_in_fsh(ifsh) +1.0e-8) ;
        phase_dlogist_fsh(ifsh)  = -1;
        logsel_dslp_in_fsh(ifsh) = 0. ;
        sel_dinf_in_fsh(ifsh)    = 0. ;        
        cout << "sel_slp_in_fsh(" << ifsh << ") = " << sel_slp_in_fsh(ifsh) << endl;
        cout << "sel_inf_in_fsh(" << ifsh << ") = " << sel_inf_in_fsh(ifsh) << endl;
      }
      break;
      case 3 : // Double logistic 
      {
        *(ad_comm::global_datafile) >> phase_sel_fsh;  
        *(ad_comm::global_datafile) >> sel_slp_in_fsh(ifsh) ;
        *(ad_comm::global_datafile) >> sel_inf_in_fsh(ifsh) ;
        *(ad_comm::global_datafile) >> sel_dslp_in_fsh(ifsh) ;
        *(ad_comm::global_datafile) >> sel_dinf_in_fsh(ifsh) ;
        for (int iyr=styr;iyr<=endyr;iyr++) { 
          *(ad_comm::global_datafile) >> sel_change_in_fsh(ifsh,iyr) ; }

        phase_selcoff_fsh(ifsh) = -1;
        phase_logist_fsh(ifsh)  = phase_sel_fsh;
        phase_dlogist_fsh(ifsh) = phase_sel_fsh+1;

        logsel_slp_in_fsh(ifsh) = log(sel_slp_in_fsh(ifsh)) ;
        logsel_dslp_in_fsh(ifsh) = log(sel_dslp_in_fsh(ifsh)) ;
        cout << "sel_slp_in_fsh(" << ifsh << ") = " << sel_slp_in_fsh(ifsh) << endl;
        cout << "sel_inf_in_fsh(" << ifsh << ") = " << sel_inf_in_fsh(ifsh) << endl;
        cout << "sel_dslp_in_fsh(" << ifsh << ") = " << sel_slp_in_fsh(ifsh) << endl;
        cout << "sel_dinf_in_fsh(" << ifsh << ") = " << sel_inf_in_fsh(ifsh) << endl;
        break;
      }
      break;
    }
  }
 END_CALCS
  ivector n_sel_ch_fsh(1,nfsh);
  imatrix yrs_sel_ch_tfsh(1,nfsh,1,endyr-styr+1);
  //!! for (ifsh=1;ifsh<=nfsh;ifsh++) 
  //!!   if (phase_selcoff_fsh(ifsh)>0) 
  //!!     curv_pen_fsh(ifsh) = 1./ (square(curv_pen_fsh(ifsh))*2);
  //!! cout << "eventual curv_pen_fsh(ifsh) = " << curv_pen_fsh(ifsh) << endl;
 LOCAL_CALCS
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
    sel_change_in_fsh(ifsh,styr)=1.;
     n_sel_ch_fsh(ifsh)=0.;
    int j=1;
    yrs_sel_ch_tfsh(ifsh,j) = styr;
    for (int iyr=styr+1;iyr<=endyr;iyr++)
    {
      if(sel_change_in_fsh(ifsh,iyr)>0)
      {
        j++;
        yrs_sel_ch_tfsh(ifsh,j) = iyr;
       cout <<"Number of select changes: fishery "<<ifsh<<": "<<n_sel_ch_fsh(ifsh)<<endl;
       cout << "yrs_sel_ch_tfsh "  << yrs_sel_ch_tfsh << endl;
      }
    }
    n_sel_ch_fsh(ifsh) = j;
    cout <<"Number of selectivity changes in fishery "<<ifsh<<": "<<n_sel_ch_fsh(ifsh)<<endl;
  }
 END_CALCS
  !!cout << "sel_change_in_fsh: "<< endl << sel_change_in_fsh << endl;
  imatrix   yrs_sel_ch_fsh(1,nfsh,1,n_sel_ch_fsh);
  imatrix   nselages_fsh(1,nfsh,1,n_sel_ch_fsh);
 LOCAL_CALCS
  for (ifsh=1; ifsh<=nfsh;ifsh++)
    yrs_sel_ch_fsh(ifsh) = yrs_sel_ch_tfsh(ifsh)(1,n_sel_ch_fsh(ifsh));
  cout<< "yrs_sel_ch_fsh: "<< endl << yrs_sel_ch_fsh << endl;
 END_CALCS

  !!cout << "bottom of DATA_SECTION" << endl << endl;
 
  // ============================
PARAMETER_SECTION
  // ============================

  !! cout << "BEGIN PARAMETER_SECTION" << endl << endl;
  //init_number        dummy(1);
  init_number dummy2(Terminal_phase);

 // Biological parameters
  init_bounded_vector log_M(1,nareas,-5,5,phase_M)
  init_bounded_matrix Mdevs(1,nareas,styr_rec,endyr,-15,15,phase_Mdevs)
  init_bounded_number log_Linf(2.5,5,phase_al);
  init_bounded_number log_vB_k(-10,5,phase_al);
  init_bounded_number log_vB_sig(0.2,5,phase_al);
  vector              wt_age(1,nages);
  number              Linf;
  number              vB_k;
  number              vB_sig;
  matrix              al_key(1,nages,1,nl_bins);
  matrix              M(1,nareas,styr_rec,endyr);
  3darray             natage(1,nareas,styr,endyr,1,nages)
  vector              Nmid(1,nages)                                             // natage mid-year for Pope's
  4darray             n_areasrv(1,nsrv,1,nareas,styr,endyr,1,nages)             // vulnerable survey biomass
  matrix              Sp_Biom(1,nareas,styr_sp,endyr)
  matrix              pred_rec(1,nareas,styr_rec,endyr)                         // from spawner-recruit relationship
  matrix              mod_rec(1,nareas,styr_rec,endyr)                          // As estimated by model
  3darray             Z(1,nareas,styr_rec,endyr,1,nages)
  3darray             S(1,nareas,styr,end_datyrs,1,nages)
  matrix              surv(1,nareas,styr_rec,endyr)
  3darray             surv_age(1,nareas,end_datyrs,endyr,1,nages)
  matrix              natmort(1,nareas,styr_rec,endyr)


 // Stock recruitment parameters
  init_vector         mean_log_rec(1,nareas,phase_meanRec); 
  init_bounded_number steepness(0.21,Steepness_UB,phase_srec)
  init_vector         log_Rzero(1,nareas,phase_Rzero)
  likeprof_number     log_Rzero_prof;         
  init_bounded_matrix rec_dev(1,nareas,styr_rec,end_datyrs,-15,15,phase_recdevs) 
  init_vector         log_sigmar(1,nareas,phase_sigmar);
  vector              alpha(1,nareas)
  vector              beta(1,nareas)
  vector              Bzero(1,nareas)                                           // pre-exploitation spawning biomass
  vector              Btot(1,nareas)                                            // prre-exploitation total biomass
  vector              Rzero(1,nareas)
  vector              phizero(1,nareas)
  vector              sigmar(1,nareas)                                          // recruitment variability
  vector              sigmarsq(1,nareas)                                        // Sigma(R)^2
 // Survey selectivity parameters
  !! for (isrv=1;isrv<=nsrv;isrv++){
  !!   msrv = msrv_id(isrv);
  !!   nselages_srv(isrv)=nselages_in_srv(msrv);  
  !! }
  init_number_vector log_q_srv(1,nsrv,phase_q_srv)                                           
  init_matrix_vector log_selcoffs_srv(1,nsrv,1,n_sel_ch_srv,1,nselages_srv,phase_selcoff_srv);
  init_vector_vector logsel_slope_srv(1,nsrv,1,n_sel_ch_srv,phase_logist_srv)       
  init_vector_vector logsel_dslope_srv(1,nsrv,1,n_sel_ch_srv,phase_dlogist_srv)       
  init_vector_vector sel50_srv(1,nsrv,1,n_sel_ch_srv,phase_logist_srv)
  init_vector_vector seld50_srv(1,nsrv,1,n_sel_ch_srv,phase_dlogist_srv)

  matrix             sel_slope_srv(1,nsrv,1,n_sel_ch_srv)
  matrix             sel_dslope_srv(1,nsrv,1,n_sel_ch_srv)
  3darray            log_sel_srv(1,nsrv,styr,endyr,1,nages)
  3darray            sel_srv(1,nsrv,styr,endyr,1,nages)
  matrix             avgsel_srv(1,nsrv,1,n_sel_ch_srv);
 // survey predictions
  3darray            pred_srv(1,nsrv,1,nareas,styr,endyr)
  !! int             lnsrv=nsrv;
  !! int             lnareas=nareas;
  !! imatrix         lnyrs_srv_comp=nyrs_srv_comp;
  !! int             lnl_bins=nl_bins;
  !! int             lnages=nages;
  4darray            eac_srv(1,lnsrv,1,lnareas,1,lnyrs_srv_comp,1,lnages)                 // Predicted survey numbers at age
  4darray            ec_srv(1,lnsrv,1,lnareas,1,lnyrs_srv_comp,1,lnl_bins)                // Predicted catch-at-length or age (depends on composition data)
  vector             q_srv(1,nsrv)                                             // Survey q

 // Fishery selectivity parameters
  !! for (ifsh=1;ifsh<=nfsh;ifsh++) nselages_fsh(ifsh)=nselages_in_fsh(ifsh);  
  init_number_vector log_q_fsh(1,nfsh,phase_q_fsh)                                           
  init_matrix_vector log_selcoffs_fsh(1,nfsh,1,n_sel_ch_fsh,1,nselages_fsh,phase_selcoff_fsh);
  init_vector_vector logsel_slope_fsh(1,nfsh,1,n_sel_ch_fsh,phase_logist_fsh)       
  init_vector_vector logsel_dslope_fsh(1,nfsh,1,n_sel_ch_fsh,phase_dlogist_fsh)       
  init_vector_vector sel50_fsh(1,nfsh,1,n_sel_ch_fsh,phase_logist_fsh)
  init_vector_vector seld50_fsh(1,nfsh,1,n_sel_ch_fsh,phase_dlogist_fsh)
  matrix             sel_slope_fsh(1,nfsh,1,n_sel_ch_fsh)
  matrix             sel_dslope_fsh(1,nfsh,1,n_sel_ch_fsh)
  3darray            log_sel_fsh(1,nfsh,styr,endyr,1,nages)
  3darray            sel_fsh(1,nfsh,styr,endyr,1,nages)
  matrix             avgsel_fsh(1,nfsh,1,n_sel_ch_fsh);

 // fishery predictions
  !! int             lnfsh=nfsh;
  //!! imatrix         lnyrs_fsh_comp=nyrs_fsh_comp;
  //!! imatrix         lnyrs_fsh = nyrs_fsh;

 // Fishing mortality parameters
  init_vector         log_avg_fmort(1,nfsh,phase_fmort)
  init_bounded_matrix fmort_dev(1,nfsh,styr,end_datyrs,-12,8,phase_fmort)
  matrix              Fmort(1,nareas,styr,end_datyrs);                        // Annual total fully-selected F mortality
  matrix              fmort_calc(1,nfsh,styr,end_datyrs);                     // F  for catch_opt == 4
  matrix              Hrate(1,nfsh,styr,endyr)                                //Harvest Rate using Pope's approximation
  3darray             catage(1,nfsh,styr,endyr,1,nages)                       // Catch-at-age
  3darray             catage_fut(1,nfsh,end_datyrs+1,endyr,1,nages)           // Catch-at-age in future
  3darray             F(1,nfsh,styr,end_datyrs,1,nages)
  3darray             eac_fsh(1,nfsh,1,nyrs_fsh_comp,1,nages)                 // Predicted survey numbers at age
  matrix              pred_catch(1,lnfsh,styr,endyr)                          // Total catch
  3darray             ec_fsh(1,nfsh,1,nyrs_fsh_comp,1,nl_bins)                // Predicted catch-at-length or age (depends on composition data)
  vector              q_fsh(1,nfsh)                                           // Survey q

 // movement among areas
  init_vector        mean_log_move(1,nmov,phase_move_m)
  init_matrix        log_move_dev(1,nmov,styr_rec,endyr,phase_move_d)
  matrix             nu(1,nmov,styr_rec,endyr)
  matrix             emig_rate(1,nareas,styr_rec,endyr)
  3darray            immigration(1,nareas,styr_rec,endyr,1,nages)
  init_bounded_vector log_influx(1,nareas,-20,5)
  vector             influx(1,nareas)

 // Likelihood and penalty value names         
  vector             rec_pen(1,4)
  number             move_pen
  vector             catch_like(1,nfsh)                                        // Catch likelihood for period with data 
  vector             comp_like_srv(1,nsrv)
  vector             comp_like_fsh(1,nfsh)
  number             CrashPen                                                  // Methot-style penalty for Pope's approximation
  matrix             selcofs_pen_srv(1,nsrv,1,4)       
  matrix             surv_like(1,nsrv,1,nareas)
  matrix             fpen(1,nareas,1,2)                                        // F-penalty
  vector             post_priors(1,4)                                          // Compute_priors penalty
  vector             post_priors_srvq(1,nsrv)                                  // Compute_priors penalty
  vector             obj_comps(1,12)
  objective_function_value obj_fun;                                            // Objective function
  !! cout << "end PARAMETER_SECTION" << endl << endl;

  sdreport_matrix SSBOut(1,nareas,styr_sp,endyr);
  sdreport_vector Bzout(1,nareas);
  sdreport_vector Btout(1,nareas);

// ==============================
INITIALIZATION_SECTION
// ==============================

  log_M             natmortprior;
  Mdevs             Mdevs_init;
  mean_log_rec      mean_log_rec_init;
  steepness         steepnessprior;
  log_Rzero         log_Rzero_init;
  logsel_slope_srv  logsel_slp_in_srv ;
  sel50_srv         sel_inf_in_srv;
  logsel_dslope_srv logsel_dslp_in_srv ;
  seld50_srv        sel_dinf_in_srv ;
  logsel_slope_fsh  logsel_slp_in_srv ;
  sel50_fsh         sel_inf_in_srv;
  logsel_dslope_fsh logsel_dslp_in_srv ;
  seld50_fsh        sel_dinf_in_srv ;
  log_influx        log_influx_init;
  log_Linf           4.1;     //60.8
  log_vB_k          -0.7985; // 0.45
  log_vB_sig        -1.897;  //0.15
  log_avg_fmort     -8.;
  fmort_dev         -8.;

// ==============================
PRELIMINARY_CALCS_SECTION
// ==============================

  cout << "Begin PRELIMINARY_CALCS_SECTION" << endl;
  log_sigmar = log(sigmarprior);
  if (DebugOut > 0 ) cout << "natmortprior = " << natmortprior << endl;
  for (isrv = 1;isrv <= nsrv;isrv++)
   for (iyr=1;iyr<=n_sel_ch_srv(isrv);iyr++)
    yrs_sel_ch_srv(isrv,iyr) = yrs_sel_ch_tsrv(isrv,iyr);
  offset_srv.initialize();
  iarea = 0;
  for(isrv=1;isrv<=nsrv;isrv++)
  for (iarea=1;iarea<=nareas;iarea++)
   {
      cout << "offset_srv(" << isrv << ") = " << endl;
      for (iyr=1;iyr<=nyrs_srv_comp(isrv,iarea);iyr++)
       {
        oc_srv(isrv,iarea,iyr) /= sum(oc_srv(isrv,iarea,iyr)); // Normalize to sum to one
        offset_srv(isrv,iarea) -= nsmpl_srv(isrv,iarea,iyr)*(oc_srv(isrv,iarea,iyr) + 0.0001) * 
                             log(oc_srv(isrv,iarea,iyr) + 0.0001 ) ;
        cout << offset_srv(isrv) << " ";
       }
      cout  << endl;
    }
  for(isrv=1;isrv<=nsrv;isrv++) log_q_srv(isrv) = log_qprior(isrv);

  for (ifsh = 1;ifsh <= nfsh;ifsh++)
   for (iyr=1;iyr<=n_sel_ch_fsh(ifsh);iyr++)
    yrs_sel_ch_fsh(ifsh,iyr) = yrs_sel_ch_tfsh(ifsh,iyr);
  offset_fsh.initialize();
  iarea = 0;
  for(ifsh=1;ifsh<=nfsh;ifsh++)
   {
      cout << "offset_fsh(" << ifsh << ") = " << endl;
      for (iyr=1;iyr<=nyrs_fsh_comp(ifsh);iyr++)
       {
        oc_fsh(ifsh,iyr) /= sum(oc_fsh(ifsh,iyr)); // Normalize to sum to one
        offset_fsh(ifsh) -= nsmpl_fsh(ifsh,iyr)*(oc_fsh(ifsh,iyr) + 0.0001) * 
                             log(oc_fsh(ifsh,iyr) + 0.0001 ) ;
        cout << offset_fsh(ifsh) << " ";
       }
      cout  << endl;
    }
  for(ifsh=1;ifsh<=nfsh;ifsh++) log_q_fsh(ifsh) = log_qprior(ifsh);
  for (imov=1;imov<=nmov;imov++)
      mean_log_move(imov) =    log_mov_init(imov);
  for (iarea=1;iarea<=nareas;iarea++)
    for(iyr=styr_rec;iyr<=end_datyrs;iyr++)
      {
        emig_rate(iarea,iyr) = 0;
        rec_dev(iarea,iyr) = rec_dev_init(iarea,iyr);  
      }
  wt_age = wt_age_init;
  cout << "End PRELIMINARY_CALCS_SECTION" << endl;
  
  // ============================
PROCEDURE_SECTION
  // ============================

  if (DebugOut > 0 ) 
    cout << endl << endl << "BEGIN PROCEDURE_SECTION" << endl << endl;
  cout << "current phase = " << current_phase() << endl;

  CrashPen = 0; // used for negative natage, catage with Pope's approximation

  log_Rzero_prof = log_Rzero(1);

  for (isrv=1;isrv<=nsrv;isrv++)
    q_srv[isrv] =  mfexp(log_q_srv[isrv]);
  for (ifsh=1;ifsh<=nfsh;ifsh++)
    q_fsh[ifsh] = mfexp(log_q_fsh[ifsh]);
  Get_Selectivity();
  Get_Mortality();
  Get_Immigration();
  Get_Bzero();
  Get_Numbers_at_Age();
  Get_Survey_Predictions();
  Catch_at_Age();
  evaluate_the_objective_function();
  for(iarea=1;iarea<=nareas;iarea++){
      Bzout(iarea) = Bzero(iarea);
      Btout(iarea) = Btot(iarea);
    for(iyr=styr_sp;iyr<=endyr;iyr++)
      SSBOut(iarea,iyr) = Sp_Biom(iarea,iyr);
      }
  if (mceval_phase())
    write_mceval();

  // ============================
FUNCTION evaluate_the_objective_function
  // ============================

  obj_comps.initialize();
  cout << endl << "begin eval. obj. function Iteration " << iters << endl;
  iters = iters + 1;
  //dvariable Temp_obj;
  //obj_fun     += dummy*dummy;

  if(catch_opt == 1)
    Cat_Like();
  Rec_Pen();
  Move_Pen();
  Age_Length();
  Comp_Like();
  Srv_Like();
  Fmort_Pen();
  Selcofs_Pen();
  if (use_priors == 1) Compute_priors();
  if (active(log_Rzero))
    for(iarea=1;iarea<=nareas;iarea++){
      obj_fun += .5 * square(log_Rzero(iarea)-mean_log_rec(iarea)); // Keep Rzero close to mean recruitment
      if (DebugOut > 0 ){
        cout << "Rzero, mean_rec SSE penalty area(" << iarea << ") = ";
        cout << .5 * square(log_Rzero(iarea)-mean_log_rec(iarea)) << endl;
        }
      }
  obj_comps(1) = sum(catch_like);
  obj_comps(2) = move_pen;
  obj_comps(3) = sum(surv_like);
  obj_comps(4) = sum(comp_like_srv);
  obj_comps(5) = sum(comp_like_fsh);
  obj_comps(6) = sum(selcofs_pen_srv);
  obj_comps(7) = sum(rec_pen);
  obj_comps(8) = sum(fpen);
  if(use_priors == 1) obj_comps(9) = sum(post_priors_srvq);
  if(use_priors == 1) obj_comps(10)= sum(post_priors);
  obj_fun     += sum(obj_comps);
  obj_fun     += dummy2*dummy2;
  obj_fun += CrashPen;

  cout << endl << "bottom eval. obj. function Iteration " << iters << endl;
  if (DebugOut > 0 ) {
    cout << endl << "mean_log_move = " << endl << mean_log_move << endl;
    cout << "3) move_pen = " << move_pen << endl;
    cout << "4a) surv_like = " << surv_like << endl;
    cout << "5) comp_like_srv = " << comp_like_srv << endl;
    cout << "6) selcofs_pen_srv = " << endl << selcofs_pen_srv << endl;
    cout << "7) rec_pen = " << rec_pen << endl;
    cout << "8) fpen = " << fpen << endl;
    cout << "obj_comps(1) (sum(catch_like) = " << obj_comps(1) << endl;
    cout << "obj_comps(2) (move_pen) = " << obj_comps(2) << endl;
    cout << "obj_comps(3) (sum(surv_like)) = " << obj_comps(3) << endl;
    cout << "obj_comps(4) (sum(comp_like_srv)) = " << obj_comps(4) << endl;
    cout << "obj_comps(5) (sum(comp_like_fsh)) = " << obj_comps(5) << endl;
    cout << "obj_comps(6) (sum(selcofs_pen_srv))= " << obj_comps(6) << endl;
    cout << "obj_comps(7) (sum(rec_pen)) = " << obj_comps(7) << endl;
    cout << "obj_comps(8) (sum(fpen)) = " << obj_comps(8) << endl;
    cout << "obj_comps(9) (sum(post_priors_srvq))= " << obj_comps(9) << endl;
    cout << "obj_comps(10) (sum(post_priors)) = " << obj_comps(10) << endl;
    cout << "CrashPen = " << CrashPen << endl;
    cout << "obj_fun = " << obj_fun << endl;
    }


  // ============================
FUNCTION Get_Selectivity
  // ============================

  if (DebugOut > 0 ) cout << "Begin Get_Selectivity" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
  {
    msrv = msrv_id(isrv);
    if (DebugOut > 0 ){
      cout << "isrv = " << isrv << endl;
      cout << "msrv = " << msrv << endl;
      }
    
    switch(srv_sel_opt(isrv))
     {
      case 1 : // Selectivity coefficients
      if (active(log_selcoffs_srv(isrv)))
      {
        int isel_ch_tmp = 1 ;
        dvar_vector sel_coffs_tmp(1,nselages_srv(isrv,isel_ch_tmp));
        for (int iyr=styr;iyr<=end_datyrs;iyr++)
         {
          if (i==yrs_sel_ch_srv(isrv,isel_ch_tmp)) 
          {
            sel_coffs_tmp.initialize();
            sel_coffs_tmp = log_selcoffs_srv(isrv,isel_ch_tmp);
            avgsel_srv(isrv,isel_ch_tmp)  = log(mean(mfexp(sel_coffs_tmp(q_age_min,q_age_max))) + 1.0e-8);
            if (isel_ch_tmp<n_sel_ch_srv(isrv))
              isel_ch_tmp++;
          }
          log_sel_srv(isrv,iyr)(1,nselages_srv(isrv,isel_ch_tmp)) = sel_coffs_tmp;
          log_sel_srv(isrv,iyr)(nselages_srv(isrv,isel_ch_tmp),nages) = 
                                 log_sel_srv(isrv,iyr,nselages_srv(isrv,isel_ch_tmp));
          log_sel_srv(isrv,iyr) -= log(mean(mfexp(log_sel_srv(isrv,iyr)(q_age_min,q_age_max))) + 1.0e-8); 
        }
       }
       break;
       case 2: // asymptotic logistic
       {
          sel_slope_srv(isrv)     = mfexp(logsel_slope_srv(msrv));
          int isel_ch_tmp         = 1 ;
          dvariable sel_slope_tmp = sel_slope_srv(isrv,isel_ch_tmp);
          dvariable sel50_tmp     = sel50_srv(msrv,isel_ch_tmp);
          for (iyr=styr;iyr<=endyr;iyr++)
          {
            if (iyr==yrs_sel_ch_srv(isrv,isel_ch_tmp)) 
            {
              sel_slope_tmp       = sel_slope_srv(isrv,isel_ch_tmp);
              sel50_tmp           = sel50_srv(isrv,isel_ch_tmp);
              if (isel_ch_tmp<n_sel_ch_srv(isrv))
                isel_ch_tmp++;
            }
            log_sel_srv(isrv,iyr)(1,nselages_srv(isrv,isel_ch_tmp)) = 
                                   -1.*log( 1.0 + mfexp(-1.*sel_slope_tmp * 
                                   (age_vector(1,nselages_srv(isrv,isel_ch_tmp)) - sel50_tmp) ) + 1.0e-8);
            log_sel_srv(isrv,iyr)(nselages_srv(isrv,isel_ch_tmp),nages) = 
                                   log_sel_srv(isrv,iyr,nselages_srv(isrv,isel_ch_tmp)) + 1.0e-8;
            log_sel_srv(isrv,iyr) -= max(log_sel_srv(isrv,iyr));  
          }
       }
       break;
      case 3 : // Double logistic
        {
          sel_slope_srv(isrv)  = mfexp(logsel_slope_srv(isrv));
          sel_dslope_srv(isrv) = mfexp(logsel_dslope_srv(isrv));
          int isel_ch_tmp = 1 ;
          dvariable sel_slope_tmp = sel_slope_srv(isrv,isel_ch_tmp);
          dvariable sel50_tmp     = sel50_srv(isrv,isel_ch_tmp);
          dvariable sel_dslope_tmp = sel_dslope_srv(isrv,isel_ch_tmp);
          dvariable seld50_tmp     = sel50_srv(isrv,isel_ch_tmp)+seld50_srv(isrv,isel_ch_tmp);
          for (iyr=styr;iyr<=endyr;iyr++)
          {
            if (iyr==yrs_sel_ch_srv(isrv,isel_ch_tmp)) 
            {
              sel_slope_tmp  = sel_slope_srv(isrv,isel_ch_tmp);
              sel50_tmp      = sel50_srv(isrv,isel_ch_tmp);
              sel_dslope_tmp = sel_dslope_srv(isrv,isel_ch_tmp);
              seld50_tmp     = sel50_srv(isrv,isel_ch_tmp)+seld50_srv(isrv,isel_ch_tmp);
              if (isel_ch_tmp<n_sel_ch_srv(isrv))
                isel_ch_tmp++;
            }
            log_sel_srv(isrv,iyr)(1,nselages_srv(isrv,isel_ch_tmp))     =
                         log( 1.0 + mfexp(-1.*sel_slope_tmp * 
                         ( age_vector(1,nselages_srv(isrv,isel_ch_tmp)) - sel50_tmp) ) + 1.0e-8)+
                         log( 1. - 1/(1.0 + mfexp(-sel_dslope_tmp * 
                         ( age_vector(1,nselages_srv(isrv,isel_ch_tmp)) - seld50_tmp))) + 1.0e-8 );

            log_sel_srv(isrv,iyr)(nselages_srv(isrv,isel_ch_tmp),nages) = 
                         log_sel_srv(isrv,iyr,nselages_srv(isrv,isel_ch_tmp));

            log_sel_srv(isrv,iyr) -= max(log_sel_srv(isrv,iyr));  
          }
        }
      break;
      case 4 : // Double normal
        {
          sel_slope_srv(isrv)  = mfexp(logsel_slope_srv(isrv));
          sel_dslope_srv(isrv) = mfexp(logsel_dslope_srv(isrv));
          //sel50_srv(isrv)      = mfexp(sel50_srv(isrv)); // case 4 only, sel50_srv already logged
          if (DebugOut > 0 ) {
            cout << "sel_slope_srv(" << isrv << ") = " << sel_slope_srv(isrv) << endl;
            cout << "sel_dslope_srv(" << isrv << ") = " << sel_dslope_srv(isrv) << endl;
            cout << "sel50_srv(" << isrv << ") = " << sel50_srv(isrv) << endl;
            }
          int isel_ch_tmp = 1 ;
          dvariable a1     = mfexp(sel50_srv(isrv,isel_ch_tmp));
          dvariable sL     = sel_slope_srv(isrv,isel_ch_tmp);
          dvariable sR     = sel_dslope_srv(isrv,isel_ch_tmp);
          for (iyr=styr;iyr<=endyr;iyr++)
          {
            if (iyr==yrs_sel_ch_srv(isrv,isel_ch_tmp)) 
            {
              a1  = mfexp(sel50_srv(isrv,isel_ch_tmp));
              sL  = sel_slope_srv(isrv,isel_ch_tmp);
              sR  = sel_dslope_srv(isrv,isel_ch_tmp);
              if (isel_ch_tmp<n_sel_ch_srv(isrv))
                isel_ch_tmp++;
            }
             for (age =1;age<nages;age++)
             log_sel_srv(isrv,iyr)(1,nselages_srv(isrv,isel_ch_tmp))(age)     =
                 log(pow(2,(-1 *pow((age_vector(1,nselages_srv(isrv,isel_ch_tmp))(age)-a1)/sL,2))) + 1.0e-8)*(1.0-
                 (-1.*log(1.0+mfexp(-1.*100*(age_vector(1,nselages_srv(isrv,isel_ch_tmp))(age)-a1))+1.0e-8))) +
                 log(pow(2,(-1 *pow((age_vector(1,nselages_srv(isrv,isel_ch_tmp))(age)-a1)/sR,2))) + 1.0e-8)*
                 (-1.*log(1.0+mfexp(-1.*100*(age_vector(1,nselages_srv(isrv,isel_ch_tmp))(age)-a1))+1.0e-8));

            log_sel_srv(isrv,iyr)(nselages_srv(isrv,isel_ch_tmp),nages) = 
                         log_sel_srv(isrv,iyr,nselages_srv(isrv,isel_ch_tmp));

            log_sel_srv(isrv,iyr) -= max(log_sel_srv(isrv,iyr));  
          }
        }
      break;
      }
  }
  sel_srv = mfexp(log_sel_srv);

  // FISHERIES SELECTIVITY
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
    switch(fsh_sel_opt(ifsh))
     {
      case 1 : // Selectivity coefficients
      if (active(log_selcoffs_fsh(ifsh)))
      {
        int isel_ch_tmp = 1 ;
        dvar_vector sel_coffs_tmp(1,nselages_fsh(ifsh,isel_ch_tmp));
        for (int iyr=styr;iyr<=end_datyrs;iyr++)
         {
          if (i==yrs_sel_ch_fsh(ifsh,isel_ch_tmp)) 
          {
            sel_coffs_tmp.initialize();
            sel_coffs_tmp = log_selcoffs_fsh(ifsh,isel_ch_tmp);
            avgsel_fsh(ifsh,isel_ch_tmp)  = log(mean(mfexp(sel_coffs_tmp(q_age_min,q_age_max))) + 1.0e-8);
            if (isel_ch_tmp<n_sel_ch_fsh(ifsh))
              isel_ch_tmp++;
          }
          log_sel_fsh(ifsh,iyr)(1,nselages_fsh(ifsh,isel_ch_tmp)) = sel_coffs_tmp;
          log_sel_fsh(ifsh,iyr)(nselages_fsh(ifsh,isel_ch_tmp),nages) = 
                                 log_sel_fsh(ifsh,iyr,nselages_fsh(ifsh,isel_ch_tmp));
          log_sel_fsh(ifsh,iyr) -= log(mean(mfexp(log_sel_fsh(ifsh,iyr)(q_age_min,q_age_max))) + 1.0e-8); 
        }
       }
       break;
       case 2: // asymptotic logistic
       {
          sel_slope_fsh(ifsh)     = mfexp(logsel_slope_fsh(ifsh));
          int isel_ch_tmp         = 1 ;
          dvariable sel_slope_tmp = sel_slope_fsh(ifsh,isel_ch_tmp);
          dvariable sel50_tmp     = sel50_fsh(ifsh,isel_ch_tmp);
           for (iyr=styr;iyr<=endyr;iyr++)
          {
            if (iyr==yrs_sel_ch_fsh(ifsh,isel_ch_tmp)) 
            {
              sel_slope_tmp       = sel_slope_fsh(ifsh,isel_ch_tmp);
              sel50_tmp           = sel50_fsh(ifsh,isel_ch_tmp);
              if (isel_ch_tmp<n_sel_ch_fsh(ifsh))
                isel_ch_tmp++;
            }
            log_sel_fsh(ifsh,iyr)(1,nselages_fsh(ifsh,isel_ch_tmp)) = 
                                   -1.*log( 1.0 + mfexp(-1.*sel_slope_tmp * 
                                   (age_vector(1,nselages_fsh(ifsh,isel_ch_tmp)) - sel50_tmp) ) + 1.0e-8);
            log_sel_fsh(ifsh,iyr)(nselages_fsh(ifsh,isel_ch_tmp),nages) = 
                                   log_sel_fsh(ifsh,iyr,nselages_fsh(ifsh,isel_ch_tmp)) + 1.0e-8;
          }
       }
       break;
      case 3 : // Double logistic
        {
          sel_slope_fsh(ifsh)  = mfexp(logsel_slope_fsh(ifsh));
          sel_dslope_fsh(ifsh) = mfexp(logsel_dslope_fsh(ifsh));
          int isel_ch_tmp = 1 ;
          dvariable sel_slope_tmp = sel_slope_fsh(ifsh,isel_ch_tmp);
          dvariable sel50_tmp     = sel50_fsh(ifsh,isel_ch_tmp);
          dvariable sel_dslope_tmp = sel_dslope_fsh(ifsh,isel_ch_tmp);
          dvariable seld50_tmp     = sel50_fsh(ifsh,isel_ch_tmp)+seld50_fsh(ifsh,isel_ch_tmp);
          for (iyr=styr;iyr<=endyr;iyr++)
          {
            if (iyr==yrs_sel_ch_fsh(ifsh,isel_ch_tmp)) 
            {
              sel_slope_tmp  = sel_slope_fsh(ifsh,isel_ch_tmp);
              sel50_tmp      =     sel50_fsh(ifsh,isel_ch_tmp);
              sel_dslope_tmp = sel_dslope_fsh(ifsh,isel_ch_tmp);
              seld50_tmp     = sel50_fsh(ifsh,isel_ch_tmp)+seld50_fsh(ifsh,isel_ch_tmp);
              if (isel_ch_tmp<n_sel_ch_fsh(ifsh))
                isel_ch_tmp++;
            }
            log_sel_fsh(ifsh,iyr)(1,nselages_fsh(ifsh,isel_ch_tmp))     =
                         -log( 1.0 + mfexp(-1.*sel_slope_tmp * 
                         ( age_vector(1,nselages_fsh(ifsh,isel_ch_tmp)) - sel50_tmp) ) + 1.0e-8)+
                         log( 1. - 1/(1.0 + mfexp(-sel_dslope_tmp * 
                         ( age_vector(1,nselages_fsh(ifsh,isel_ch_tmp)) - seld50_tmp))) + 1.0e-8 );

            log_sel_fsh(ifsh,iyr)(nselages_fsh(ifsh,isel_ch_tmp),nages) = 
                         log_sel_fsh(ifsh,iyr,nselages_fsh(ifsh,isel_ch_tmp));

              log_sel_fsh(ifsh,iyr) -= max(log_sel_fsh(ifsh,iyr));  
          }
        }
      break;
      }
  }
  sel_fsh = mfexp(log_sel_fsh);
  if (DebugOut > 0 ) cout << "End Get_Selectivity" << endl;
  
  // ============================
FUNCTION Get_Mortality
  // ============================

  if (DebugOut > 0 ) cout << "begin Get_Mortality" << endl;
  surv.initialize();
  S.initialize();
  nu.initialize();
  emig_rate.initialize();
  dvariable Ctmp; // Baranov for catch_opt == 4
  dvariable Dtmp; // (Baranov - pred_catch) differentiated wrt F for catch_opt == 4
  //F.initialize();
  imov = 0;
  for(iarea=1;iarea<=nareas;iarea++)
    for(iyr=styr_rec;iyr<=endyr;iyr++)
      M(iarea,iyr) = mfexp(log_M(iarea)); // + mfexp(Mdevs(iarea,iyr));
    for(iarea=1;iarea<=nareas;iarea++)
      for(jarea=1;jarea<=(nareas-1);jarea++)
       {
        imov = imov + 1;        
        for(iyr=styr_rec;iyr<=endyr;iyr++)
         { 
          nu(imov,iyr) = mean_log_move(imov) + log_move_dev(imov,iyr);
          emig_rate(iarea,iyr) += mfexp(nu(imov,iyr));
         } 
       }
  for(iarea=1;iarea<=nareas;iarea++){
    for(iyr=styr_rec;iyr<=endyr;iyr++){
      natmort(iarea,iyr) = M(iarea,iyr)+emig_rate(iarea,iyr);
      Z(iarea,iyr)       = M(iarea,iyr)+emig_rate(iarea,iyr);
      surv(iarea,iyr)  = mfexp(-1.0* (M(iarea,iyr)+emig_rate(iarea,iyr)))+1e-8;
      }
    for(iyr=end_datyrs;iyr<=endyr;iyr++)
      surv_age(iarea,iyr)     = mfexp(-1*M(iarea,iyr)+emig_rate(iarea,iyr))+1e-8;
    }   // end iarea
  Fmort.initialize();
    for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      iarea = areas_fsh(ifsh);
      if (catch_opt == 1 || catch_opt == 2) // Ianelli F or catches known
        {
         Fmort(iarea) += mfexp(log_avg_fmort(ifsh) +  fmort_dev(ifsh))+1.0e-8;
          for (iyr=styr;iyr<=end_datyrs;iyr++)
            F(ifsh,iyr)   = mfexp(log_avg_fmort(ifsh) + fmort_dev(ifsh,iyr)) * sel_fsh(ifsh,iyr)+1.0e-8;
        }

      if (catch_opt == 4)       // Martell method
         {
        // CrashPen.initialize();
	 dvar_matrix Ztmp(styr,end_datyrs,1,nages);
          if (DebugOut > 0) cout << "catch_opt == 4" << endl;
          fmort_calc.initialize();
          for (iyr=styr;iyr<end_datyrs;iyr++)
           {
            fmort_calc(ifsh,iyr)   = obs_catch(ifsh,iyr) / (mfexp(-1*natmort(iarea,iyr)/2) *
                        sum(elem_prod(elem_prod(sel_fsh(ifsh,iyr),wt_age),natage(iarea,iyr))+1.0e-8));
	    if (DebugOut > 0 ) cout << "fmort_calc(ifsh) = " << fmort_calc << endl;
           for(i=1;i<6;i++) // iterate 5 times to match fmort_calc to obs_catch
            {
            Ztmp(iyr)    = fmort_calc(ifsh,iyr)*sel_fsh(ifsh,iyr) + natmort(iarea,iyr);

             Ctmp           = sum(elem_prod(elem_div(elem_prod(elem_prod(sel_fsh(ifsh,iyr),wt_age),
                                natage(iarea,iyr))*fmort_calc(ifsh,iyr),Ztmp(iyr)+1.0e-8),(1-mfexp(-Ztmp(iyr))))
                               );

             Dtmp           = sum( elem_div(elem_prod(elem_prod(natage(iarea,iyr),wt_age),
                               (natmort(iarea,iyr)+elem_prod(mfexp(-Ztmp(iyr)),
                               (elem_prod(fmort_calc(ifsh,iyr)*sel_fsh(ifsh,iyr),
                                Ztmp(iyr))-natmort(iarea,iyr))))),elem_prod(Ztmp(iyr),Ztmp(iyr)))
                               );
            fmort_calc(ifsh,iyr) = fmort_calc(ifsh,iyr) - (Ctmp-obs_catch(ifsh,iyr)) / (Dtmp+1.0e-8);
            } // end i loop
            if (DebugOut > 0 ) {
              cout << "Ctmp(" << iyr << ") = " << Ctmp << endl;
              cout << "obs_catch(" << iyr << ") = " << obs_catch(ifsh,iyr) << endl;
              }
           CrashPen += sqrt((Ctmp-obs_catch(ifsh,iyr))*(Ctmp-obs_catch(ifsh,iyr))+1.0e-8);
           } // end iyr
          }  // end catch_opt == 4 
	if(catch_opt==3)       // Pope CrashPen when catage > natage
        {
         for (iyr=styr;iyr<end_datyrs;iyr++)
         {
          ifsh = iarea;
           for (age=nages;age>1;age--) // if not enough krill of age to supply catch, Methot-style penalty
            {
             ifsh = iarea;
             if(catage(ifsh,iyr,age) > natage(ifsh,iyr,age)) // not enough krill of age to supply catch
               CrashPen += 100*catage(ifsh,iyr,age)/natage(ifsh,iyr,age);
            }
          } // end iyrs
         }  // end CrashPen
        for (iyr=styr;iyr<=end_datyrs;iyr++)
         {
          if (catch_opt != 4) Z(iarea,iyr)     += F(ifsh,iyr);
          else                Z(iarea,iyr)     += fmort_calc(ifsh,iyr);
          S(iarea,iyr)     = mfexp(-1*Z(iarea,iyr))+1e-8;
         }
       } // end ifsh
       if (DebugOut > 0 ) cout << "F = " << endl << F << endl;

  if (DebugOut > 0 ) cout << "End Get_Mortality" << endl;

  // ============================
FUNCTION Get_Immigration
  // ============================

  if (DebugOut > 0 ) cout << "begin Get_Immigration" << endl;
  immigration.initialize();
  imov = 0; 
  for (iarea=1;iarea<=nareas;iarea++)
   {
    for(jarea=1;jarea<=nareas;jarea++)
    {
      // styr_rec to styr
      if (jarea != iarea) 
      {
      imov = imov + 1;
      for(iyr=styr_rec;iyr < styr;iyr++)
        {
         immigration(iarea,iyr) += elem_prod(
            elem_prod(natage(jarea,styr),(1.0 - mfexp(-Z(jarea,iyr)))),
            mfexp(nu(imov,iyr))/(Z(jarea,iyr)));
        }  
      }
      for(iyr=styr;iyr<=endyr;iyr++)
        if (jarea != iarea)
        {
        // styr to endyr
        immigration(iarea,iyr) += elem_prod(
            elem_prod(natage(jarea,iyr),(1.0 - mfexp(-Z(jarea,iyr)))),
            mfexp(nu(imov,iyr))/Z(jarea,iyr));         
        }
        if (current_phase() < phase_move_m) immigration(iarea) = 0;
     }
    }
  if (DebugOut > 0 ) cout << "End Get_Immigration" << endl;

  // ============================
FUNCTION Get_Bzero
  // ============================

  if (DebugOut > 0 ) cout << "Begin Get_Bzero" << endl;
  Bzero.initialize();
  Sp_Biom.initialize();
  if (influx_on > 0)
      influx.initialize();
  else log_influx = 0;
  dvariable   survtmp;
  dvar_matrix natagetmp(styr_rec,styr,1,nages);
  iarea = 0;
  for (iarea=1;iarea<=nareas;iarea++)
  {
      natagetmp.initialize();
      survtmp               = mfexp(-(natmort(iarea,styr)));
      if (DebugOut > 0 ) cout << "survtmp = " << survtmp << endl;
      Rzero(iarea)          = mfexp(log_Rzero(iarea)); 
      natagetmp(styr_rec,1) = Rzero(iarea);
    if (influx_on > 0)
      influx(iarea)         = mfexp(log_influx(iarea));
    if (influx_on == 0)
      influx(iarea)         = 0;
      for (j=2; j<=nages; j++)
       {
        natagetmp(styr_rec,j) = ((natagetmp(styr_rec,j-1) + 
                                immigration(iarea,styr,j-1)) * survtmp)+
                                ((natagetmp(styr_rec,j-1) + 
                                immigration(iarea,styr,j-1)) * influx(iarea) * survtmp);
       }
      natagetmp(styr_rec,nages) /= (1.-survtmp +1e-8);
 
      Bzero(iarea)          = elem_prod(wt_age , mat_age) * survtmp * natagetmp(styr_rec);
      Btot(iarea)           = wt_age * survtmp * natagetmp(styr_rec);
      if (DebugOut > 0 ) {
        cout << endl << "Bzero(" << iarea << ")= " << Bzero(iarea) << endl;
        cout << endl << "survtmp = " << survtmp << endl;
        cout << endl << "natagetmp(styr_rec) = " << natagetmp(styr_rec) << endl;
        cout << endl << "elem_prod(wt_age , mat_age) = " << elem_prod(wt_age , mat_age) << endl;
        }
      phizero(iarea)        = Bzero(iarea)/Rzero(iarea);
      switch (SrType)
      {
       case 1:
        alpha = log(-4.*steepness/(steepness-1.));
        break;
       case 2:
        alpha(iarea)        =  Bzero(iarea) * (1. - (steepness - 0.2) / (0.8*steepness) ) / Rzero(iarea);
        beta(iarea)         = (5. * steepness - 1.) / (4. * steepness * Rzero(iarea));
        break;
       case 4:
        //R = S * EXP(alpha - beta * S))
        beta(iarea)  = log(5.*steepness)/(0.8*Bzero(iarea)) ;
        alpha(iarea) = log(Rzero(iarea)/Bzero(iarea))+beta(iarea)*Bzero(iarea);
        break;
      }
      Sp_Biom(iarea)(styr_sp,styr_rec) = Bzero(iarea);
      for (iyr=styr_rec;iyr<styr;iyr++)
       {
        Sp_Biom(iarea,iyr)  = natagetmp(iyr) * 
                               surv(iarea,iyr) * elem_prod(wt_age,mat_age); // * influx(iarea)
        natagetmp(iyr,1)    = mfexp(rec_dev(iarea,iyr) + mean_log_rec(iarea));
        natagetmp(iyr+1)(2,nages) = ++ (natagetmp(iyr)(1,nages-1)* 
                                        mfexp(-(natmort(iarea,iyr)))); 
        natagetmp(iyr+1,nages) += natagetmp(iyr,nages)*
                                   mfexp(-(natmort(iarea,iyr)));
       }
      natagetmp(styr,1)    = mfexp(rec_dev(iarea,styr) + mean_log_rec(iarea));
      mod_rec(iarea)(styr_rec,styr) = column(natagetmp,1);
      natage(iarea,styr)    = natagetmp(styr); 
      for(isrv=1;isrv<=nsrv;isrv++)
        n_areasrv(isrv,iarea,styr) = natage(iarea,styr)*pow(S(iarea,styr,1),rho(isrv));
    }
  if (DebugOut > 0 ) {
    cout << "alpha = " << alpha << endl;
    cout << "beta = " << beta << endl;
    cout << "End Get_Bzero" << endl;
    }
 
  //=============================
FUNCTION Get_Numbers_at_Age
  //=============================

  if (DebugOut > 0 ) {
    cout << "Begin Get_Numbers_at_Age" << endl;
    cout << "mean_log_rec = " << mean_log_rec << endl;
    cout << "rec_dev = " << endl << rec_dev << endl;
    }

  iarea = 0;
  for (iarea=1;iarea<=nareas;iarea++)
     {
      ifsh = iarea;
       // Recruitment in start year
         natage(iarea,styr,1) = mfexp(mean_log_rec(iarea) + rec_dev(iarea,styr));
       // Recruitment in subsequent years
         for (iyr=styr+1;iyr<=end_datyrs;iyr++)
          natage(iarea,iyr,1)= mfexp(mean_log_rec(iarea) + rec_dev(iarea,iyr));
      mod_rec(iarea,styr)  = natage(iarea,styr,1);

      // cohort numbers at age during the data period
        if (catch_opt != 3) // catch removed continuously
        for (iyr=styr;iyr<end_datyrs;iyr++)
         {
          natage(iarea,iyr+1)(2,nages)= elem_prod(natage(iarea,iyr)(1,nages-1),
                                         S(iarea,iyr)(1,nages-1))++;
          natage(iarea,iyr+1,nages)+= natage(iarea,iyr,nages)*
                                         S(iarea,iyr,nages);
          Sp_Biom(iarea,iyr)  = elem_prod(natage(iarea,iyr),S(iarea,iyr)) * 
                                         elem_prod(wt_age,mat_age) ;
          mod_rec(iarea,iyr+1)  = natage(iarea,iyr+1,1);

          for (age=1;age<=nages;age++)
          if(natage(iarea,iyr+1,age)<=0.0)                          //   Methot-style penalty
            {
             CrashPen += log(1.0 - natage(iarea,iyr+1,age));
             natage(iarea,iyr+1,age) = natage(iarea,iyr+1,age)*10.;
            }
         }
	 
        else      
        // catch removed at midpoint of year (Pope's approximation)
        for (iyr=styr;iyr<end_datyrs;iyr++)
        {
          for(age=2;age<=nages;age++) //  following Methot SS2 tpl
            {
              if (catage(ifsh,iyr)(age-1)/(natage(iarea,iyr)(age-1)*mfexp(-0.5*M(iarea,iyr))) <= 0.95)

                natage(iarea,iyr+1)(age)= ((natage(iarea,iyr)(age-1)*mfexp(M(iarea,iyr))))- 
                                           ((catage(ifsh,iyr)(age-1))*mfexp(-0.5*M(iarea,iyr)));
             
              else
               CrashPen += (catage(ifsh,iyr)(age-1)/
                           (natage(iarea,iyr)(age-1)*mfexp(-0.5*M(iarea,iyr)))) *100.;

              if(natage(iarea,iyr+1,age)<=0.0)                          //   Methot-style penalty
               {
                CrashPen += log(1.0 - natage(iarea,iyr+1,age));
                natage(iarea,iyr+1,age) = natage(iarea,iyr+1,age)*10.;
               }
	       
           natage(iarea,iyr+1)(nages) += ((natage(iarea,iyr)(nages-1)*mfexp(M(iarea,iyr))))- 
                                          ((catage(ifsh,iyr)(nages-1))*mfexp(-0.5*M(iarea,iyr)));

           if(natage(iarea,iyr+1,nages)<=0.0)                          
             {
              CrashPen += log(1.0 - natage(iarea,iyr+1,nages));
              natage(iarea,iyr+1,nages) += ((natage(iarea,iyr)(nages-1)*mfexp(M(iarea,iyr)))+ 
                                          (catage(ifsh,iyr)(nages-1)))*mfexp(-0.5*M(iarea,iyr));;
             }
           } // end of age loop

        } // end of Pope's approximation
      } // end of iarea


  for (iarea=1;iarea<=nareas;iarea++){
        for (iyr=styr;iyr<end_datyrs;iyr++)
        {
          Sp_Biom(iarea,iyr)  = elem_prod(natage(iarea,iyr),S(iarea,iyr)) * 
                                         elem_prod(wt_age,mat_age) ;
          mod_rec(iarea,iyr+1)  = natage(iarea,iyr+1,1);
        }
      if (DebugOut > 0 ) {
        cout << "natage(iarea,end_datyrs) = " << natage(iarea,end_datyrs) << endl;
        cout << "S(iarea,end_datyrs) = " << S(iarea,end_datyrs) << endl;
        }
      Sp_Biom(iarea,end_datyrs)   = elem_prod(natage(iarea,end_datyrs),S(iarea,end_datyrs)) * 
                                         elem_prod(wt_age,mat_age);
      if (DebugOut > 0 ) {
        cout << "Sp_Biom(iarea,end_datyrs) = " << Sp_Biom(iarea,end_datyrs) << endl;
        cout << "natage - catage =" << endl << natage - catage << endl;
        }

      /////////////////////  FUTURE PROJECTIONS  /////////////////////
      if (endyr > end_datyrs)
      {                           // styr for projection period is end_datyrs+1
       ifsh = iarea;
       dvar_vector rec_dev_future(end_datyrs+1, endyr);
       rec_dev_future.initialize();
       rec_dev_future = rec_dev(iarea)(end_datyrs-(nproj_yrs+2),end_datyrs-3).shift(end_datyrs+1)(end_datyrs+1, endyr);
      // Recruitment in first projection year
         natage(iarea,end_datyrs+1,1) = ((mfexp(mean_log_rec(iarea) + rec_dev_future(end_datyrs+1)) +
                              immigration(iarea,end_datyrs+1,1))) +
                             ((mfexp(mean_log_rec(iarea) + rec_dev_future(end_datyrs+1)) +
                              immigration(iarea,end_datyrs+1,1)) * influx(iarea)); 
      // Recruitment (age 1 column) in subsequent projection years
         for (iyr=end_datyrs+2;iyr<=endyr;iyr++)
          natage(iarea,iyr,1)= (((mfexp(mean_log_rec(iarea) + rec_dev_future(iyr)) +
                               immigration(iarea,iyr,1))) + 
                               ((mfexp(mean_log_rec(iarea) + rec_dev_future(iyr)) +
                               immigration(iarea,iyr,1)) * influx(iarea)));
         mod_rec(iarea,end_datyrs+1)  = natage(iarea,end_datyrs+1,1);

      // cohort numbers at age (columns 2 to 10) during the projection period
         // first projection year
         natage(iarea,end_datyrs+1)(2,nages)= (((elem_prod((natage(iarea,end_datyrs)(1,nages-1)+
                                         immigration(iarea,end_datyrs)(1,nages-1)),
                                         surv_age(iarea,end_datyrs)(1,nages-1)))++) +
                                        (((elem_prod((natage(iarea,end_datyrs)(1,nages-1)+
                                         immigration(iarea,end_datyrs)(1,nages-1)),
                                         surv_age(iarea,end_datyrs)(1,nages-1)))++ * influx(iarea))));
        // subsequent projection years 
        for (iyr=end_datyrs+1;iyr<endyr;iyr++)        
        {
          natage(iarea,iyr+1)(2,nages)= (((elem_prod((natage(iarea,iyr)(1,nages-1)+
                                         immigration(iarea,iyr)(1,nages-1)),
                                         surv_age(iarea,iyr)(1,nages-1)))++) +
                                        (((elem_prod((natage(iarea,iyr)(1,nages-1)+
                                         immigration(iarea,iyr)(1,nages-1)),
                                         surv_age(iarea,iyr)(1,nages-1)))++ * influx(iarea))))- 
                                         catage_fut(ifsh,iyr)(1,nages-1)++;
          natage(iarea,iyr+1,nages)+= (((natage(iarea,iyr,nages))*
                                         surv_age(iarea,iyr,nages)) +
                                         ((natage(iarea,iyr,nages))*
                                         surv_age(iarea,iyr,nages) * influx(iarea))) - 
                                         catage_fut(ifsh,iyr,nages);          
          Sp_Biom(iarea,iyr)  = elem_prod(natage(iarea,iyr),surv_age(iarea,iyr)) * 
                                         elem_prod(wt_age,mat_age) ;
          mod_rec(iarea,iyr+1)  = natage(iarea,iyr+1,1);
        }
        for (iyr=end_datyrs+1;iyr<=endyr;iyr++)
          for(age=1;age<=nages;age++)                      // correct for future catage greater than natage
            if(natage(ifsh,iyr,age) < 0) natage(ifsh,iyr,age) = 0;

      } // end if (endyr > end_datyrs) //future projections
      Sp_Biom(iarea,endyr)  = natage(iarea,endyr) * elem_prod(surv_age(iarea,endyr), 
                                         elem_prod(wt_age,mat_age)); 
      if (DebugOut > 0 ) cout << "Sp_Biom(iarea,endyr) = "  << Sp_Biom(iarea,endyr) << endl;
      for(isrv=1;isrv<=nsrv;isrv++)
       {
        n_areasrv(isrv,iarea,styr,1) = natage(iarea,styr,1)*pow(S(iarea,styr,1),rho(iarea));
        for (iyr=styr+1;iyr<=end_datyrs;iyr++)
          n_areasrv(isrv,iarea,iyr,1) = natage(iarea,iyr,1)*pow(S(iarea,iyr,1),rho(iarea));
        for (iyr=styr;iyr<=end_datyrs;iyr++)
         {
          n_areasrv(isrv,iarea,iyr)(2,nages-1)=elem_prod(natage(iarea,iyr)(2,nages-1),
                                                 pow(S(iarea,iyr)(2,nages-1),rho(iarea)));
          n_areasrv(isrv,iarea,iyr,nages)=natage(iarea,iyr,nages)*
                                                 pow(S(iarea,iyr,nages),rho(iarea));
         }
         //cout << "n_areasrv(" << isrv << ") = " << endl << n_areasrv(isrv) << endl;
       }
     } // end of iarea for future projections

     if (DebugOut > 0 ) cout << "End Get_Numbers_at_Age" << endl;

  // ============================
FUNCTION Get_Survey_Predictions
  // ============================

  if (DebugOut > 0 ) cout << "Begin Get_Survey_Predictions" << endl;
  dvariable sum_tmp;
  sum_tmp.initialize();
  int yy;
  for (isrv=1;isrv<=nsrv;isrv++)
  {
   iarea = 0;
   for (iarea=1;iarea<=nareas;iarea++)
    {
      for (iyr=styr;iyr<=endyr;iyr++)
       {
          pred_srv(isrv,iarea,iyr) = q_srv(isrv) * n_areasrv(isrv,iarea,iyr) * 
                                elem_prod(sel_srv(isrv)(iyr) , wt_age);       
       }
      for (iyr=1;iyr<=nyrs_srv_comp(isrv,iarea);iyr++)
       {
        yy = yrs_srv_comp(isrv,iarea,iyr); 
        dvar_vector tmp_n = elem_prod(sel_srv(isrv,yy), n_areasrv(isrv,iarea,yy));  
        sum_tmp           = sum(tmp_n)+1.0e-10;
        eac_srv(isrv,iarea,iyr)      = tmp_n/sum_tmp;
       }
    }
  }
  if (DebugOut > 0 ){
    //cout << "eac_srv" << endl << eac_srv << endl;
    cout << "End Get_Survey_Predictions" << endl;
    }

  
  //=============================
FUNCTION Catch_at_Age
  //=============================

  if (DebugOut == 1) cout << "Begin Catch_at_Age" << endl;
  dvariable fut_catch_tmp; // future catches will be gamma * Btot(iarea);
  fut_catch_tmp.initialize();
  dvariable cat_wt_missing; // not enough fish at age to supply the catch given selectivity
  for (ifsh=1; ifsh<=nfsh; ifsh++){
      iarea = areas_fsh(ifsh);
      if(endyr > end_datyrs)
      {
        fut_catch_tmp = gamma * Btot(iarea);
        for (iyr=end_datyrs+1; iyr<=endyr; iyr++) // catches for projection period as assigned
          obs_catch(ifsh,iyr) = value(fut_catch_tmp);
      }

      for (iyr=styr;iyr<=end_datyrs;iyr++) // catage during DATA YEARS
        {
        if (catch_opt == 1) // if catch is an estimated quantity
          catage(ifsh,iyr) = q_fsh(ifsh) * elem_prod(elem_div(F(ifsh,iyr),Z(iarea,iyr)),elem_prod(1.-S(iarea,iyr),natage(iarea,iyr)));
        if (catch_opt > 1) // catch is considered known without error
          catage(ifsh,iyr) = obs_catch(ifsh,iyr) *  elem_div((elem_prod(sel_fsh(ifsh,iyr),wt_age))/
                           (1.0e-10+sum(elem_prod(sel_fsh(ifsh,iyr),wt_age))),wt_age);
        }
        //if (DebugOut > 0) cout << "catage = " << endl << catage << endl;
       //elem_prod(elem_div(F(ifsh,iyr),Z(iarea,iyr)),elem_prod(1.-S(iarea,iyr),natage(iarea,iyr)));

      if(endyr > end_datyrs)
      for (iyr=end_datyrs+1;iyr<=endyr;iyr++)  // catage during PROJECTION YEARS
        {
        catage(ifsh,iyr) = obs_catch(ifsh,iyr) * elem_div((elem_prod(sel_fsh(ifsh,iyr),wt_age))/
                           (1.0e-10+sum(elem_prod(sel_fsh(ifsh,iyr),wt_age))),wt_age);
        catage_fut(ifsh,iyr) = catage(ifsh,iyr);
        for(age=nages;age>1;age--)
          if(catage_fut(ifsh,iyr,age) > natage(ifsh,iyr,age)) // not enough krill of age to supply catch
            {
             // move future catch from older ages to younger
             // ages abundant enough to support it
             cat_wt_missing = (catage_fut(ifsh,iyr,age)-natage(ifsh,iyr,age))*wt_age(age);
             catage_fut(ifsh,iyr,age) = natage(ifsh,iyr,age);   
             catage_fut(ifsh,iyr,age-1) = catage_fut(ifsh,iyr,age-1)+
                                          cat_wt_missing/wt_age(age-1); 
             cat_wt_missing = 0;
            }       
         }


      dvar_matrix Ctmp = catage(ifsh); // Copy 3darray to matrix for efficiency...
      for (iyr=styr; iyr<=endyr; iyr++)
         pred_catch(ifsh,iyr) = Ctmp(iyr)*wt_age;
      for (iyr=1; iyr<=nyrs_fsh_comp(ifsh); iyr++)
         eac_fsh(ifsh,iyr)=Ctmp(yrs_fsh_comp(ifsh,iyr))/(1.0e-10+sum(Ctmp(yrs_fsh_comp(ifsh,iyr))));
      }  
  if (DebugOut == 1){
    cout << "eac_fsh" << endl << eac_fsh << endl;
    cout << "End Catch_at_Age" << endl;
    }
 
  // ============================
FUNCTION dvar_vector SRecruit(const dvar_vector& Stmp)
  // ============================

  RETURN_ARRAYS_INCREMENT();
  int i_area = iarea;
  dvar_vector RecTmp(Stmp.indexmin(),Stmp.indexmax());
    switch (SrType)
    {
    case 1:
      RecTmp = elem_prod((Stmp / phizero(i_area)) , mfexp( alpha(i_area) * ( 1. - Stmp / Bzero(i_area) ))) ; //Ricker form from Dorn
      break;
    case 2:
      RecTmp = elem_prod(Stmp , 1. / ( alpha(i_area) + beta(i_area) * Stmp));        //Beverton-Holt form
      break;
    case 3:
      RecTmp = Rzero(i_area)*mfexp(mean_log_rec(iarea));                             //Avg recruitment
      break;
    case 4:
      RecTmp = elem_prod(Stmp , mfexp( alpha(i_area)  - Stmp * beta(i_area))) ;      //Old Ricker form
      break;
    }  
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;

  // ============================
FUNCTION Rec_Pen
  // ============================

  if (DebugOut > 0 ) cout << "Begin Rec_Pen" << endl;  
  rec_pen.initialize(); 
  //dvar sigmar_fut; sigmar_fut.initialize();
  if (active(rec_dev))
  {
   iarea = 0;
   for (iarea=1;iarea<=nareas;iarea++)
   {
    sigmar(iarea)     =  mfexp(log_sigmar(iarea));
    sigmarsq(iarea)   =  square(sigmar(iarea));
      dvariable SSQRec;
      SSQRec.initialize();
      //if (current_phase()==phase_srec) //cout << "here0" << endl;
      if (current_phase()>phase_srec) // if steepness estimated
        {
         pred_rec(iarea) = SRecruit(Sp_Biom(iarea)(styr_rec-rec_age,endyr-rec_age).shift(styr_rec)(styr_rec,endyr));
         if (DebugOut > 0 ){
           cout << "alpha(" << iarea << ") = " << alpha(iarea) << endl;
           cout << "beta(" << iarea << ") = " << beta(iarea) << endl;
           cout << "pred_rec(" << iarea << ") = " << SRecruit(Sp_Biom(iarea)(styr_rec-rec_age,endyr-rec_age).shift(styr_rec)(styr_rec,endyr)) << endl;
           }
         dvar_vector chi = log(elem_div(mod_rec(iarea)(styr_rec_est,endyr_rec_est) ,
                                  pred_rec(iarea)(styr_rec_est,endyr_rec_est) + 1.0e-8) + 1.0e-8);
         SSQRec   = norm2( chi ) ;
         rec_pen(1) += norm2(chi + sigmarsq(iarea)/2.)/(2*sigmarsq(iarea)) + nrecs_est*log_sigmar(iarea);
        }
      if (last_phase())
       {
         // Variance term for the parts not estimated by sr curve
         rec_pen(4) += .5*norm2( rec_dev(iarea)(styr_rec,styr_rec_est) )/sigmarsq(iarea) + 
                         (styr_rec_est-styr_rec)*log(sigmar(iarea) + 1.0e-8) ; 
         //if (endyr>endyr_rec_est)
         // rec_pen(4) += .5*norm2(rec_dev(iarea)(endyr_rec_est,end_datyrs) )/sigmarsq(iarea) + 
         //                  (end_datyrs-endyr_rec_est)*log(sigmar(iarea) + 1.0e-8) ; 
       }
      else
        rec_pen(2) += norm2( rec_dev(iarea)(styr_rec_est,endyr_rec_est) ) ;
      rec_pen(2) += norm2( rec_dev(iarea)(styr_rec_est,end_datyrs) ) ;
      //for(iyr = end_datyrs;iyr <= endyr; iyr++;)
      //  rec_like(3) += norm2(rec_dev(iarea,iyr))/(2*square(sigmar_fut(iarea)))+ nproj_yrs*log(sigmar_fut);
   }
  }
     //rec_pen(3) = 10*mean(rec_dev(iarea)(styr_rec,styr-1)) - mean(rec_dev(iarea)(styr,endyr));
     if (DebugOut > 0 ) cout << "pred_rec = " << pred_rec << endl;
     if (DebugOut > 0 ) cout << "mod_rec = " << mod_rec << endl;
     if (DebugOut > 0 ) cout << "End Rec_like" << endl;  

  // ============================
FUNCTION Age_Length
  // ============================

  #include "statsLib.h"
  // al_key(1,nages,1,nl_bins)
  al_key.initialize();
  wt_age.initialize();
  Linf = mfexp(log_Linf);
  vB_k = mfexp(log_vB_k);
  vB_sig = mfexp(log_vB_sig);
  dvariable z1;
  dvariable z1a;
  dvariable z2;
  dvariable z2a;
  if (DebugOut > 0 ){
    cout << "Linf = " << Linf << endl;
    cout << "vB_k = " << vB_k << endl;
    cout << "vB_sig = " << vB_sig << endl;
    }
  int slen,nlen;slen=lbins.indexmin();nlen=lbins.indexmax()-1;
  dvar_vector al_key_m(1,nages);
  if (Linf * (1-mfexp(-vB_k)) > 10) // average length of 1 year olds must be at least 10 mm
  for (age=1;age<=nages;age++){
    al_key_m(age) = Linf * (1-mfexp(-vB_k * age));
    for (ilen=slen;ilen<(nlen+1);ilen++){
      z1 = (lbins(ilen)-al_key_m(age))/vB_sig;
      z2 = (lbins(ilen+1)-al_key_m(age))/vB_sig;
      al_key(age,ilen) = cumd_norm(z2)-cumd_norm(z1);
      }
   z1a = (lbins(nlen)-al_key_m(age))/vB_sig;
   z2a = (lbins(lbins.indexmax())-al_key_m(age))/vB_sig;
   al_key(age,lbins.indexmax()) = cumd_norm(z2a)-cumd_norm(z1a);
   al_key(age) /= sum(al_key(age));
   }  
  for (age=1;age<=nages;age++)
   for(ilen=slen;ilen<=nlen;ilen++){
     //if(al_key(age,ilen) > 1e-3)
     wt_age(age) += wt_len(ilen)*al_key(age,ilen)+1e-10;
     }


  // ============================
FUNCTION Move_Pen
  // ============================

  move_pen.initialize();
  if (active(log_move_dev))
      move_pen += norm2(log_influx ) ;

  // ============================
FUNCTION Cat_Like
  // ============================
  if (DebugOut > 0) cout << "begin Cat_Like" << endl;
  
  catch_like.initialize();
  for (ifsh=1; ifsh<=nfsh; ifsh++){
    iarea = areas_fsh(ifsh);
    if (DebugOut > 0) {
      cout << "obs_catch(ifsh) = " << endl << obs_catch(ifsh) << endl;
      cout << "obs_catch_pen(ifsh) = " << endl << obs_catch_pen(ifsh) << endl;
      cout << "pred_catch(ifsh) = " << endl << pred_catch(ifsh) << endl;
      }
    if(catch_opt == 1)
    catch_like(ifsh) = obs_catch_pen(iarea) * 
                       norm2(log(elem_prod(effort(ifsh)(styr,end_datyrs),obs_catch(ifsh)(styr,end_datyrs)) +.000001)-
                             log(elem_prod(effort(ifsh)(styr,end_datyrs),pred_catch(ifsh)(styr,end_datyrs)) +.000001)+0.01);
   }
  if (DebugOut > 0) cout << "end Cat_Like" << endl;

  // ============================
FUNCTION Comp_Like
  // ============================

  if (DebugOut > 0 ) {
    cout << "Begin Comp_Like" << endl;
    cout << "eac_srv = " << eac_srv << endl;
    }
  comp_like_srv.initialize();
  comp_like_fsh.initialize();
  for(isrv=1;isrv<=nsrv;isrv++)
   {
    for (iarea=1;iarea<=nareas;iarea++)
     if (nyrs_srv_comp(isrv,iarea) >0)
     {
      ec_srv(isrv,iarea) = eac_srv(isrv,iarea)*al_key;
      for (iyr=1;iyr<=nyrs_srv_comp(isrv,iarea);iyr++){
       comp_like_srv(isrv) -= nsmpl_srv(isrv,iarea,iyr)*(oc_srv(isrv,iarea,iyr) + 0.0001) * 
                              log(ec_srv(isrv,iarea,iyr)  + 0.0001 ) ;
       }
     }
     if (DebugOut > 0 ) {
       cout << "comp_like_srv(" << isrv << ") before offset:" << endl;
       cout << comp_like_srv << endl;
       comp_like_srv(isrv) -= offset_srv(isrv);
       cout << "comp_like_srv(" << isrv <<") after offset:" << endl;
       cout << comp_like_srv << endl;
       }
    }

  for(ifsh=1;ifsh<=nfsh;ifsh++)
   {
     if (nyrs_fsh_comp(ifsh) >0)
     {
      iarea = areas_fsh(ifsh);
      ec_fsh(ifsh) = eac_fsh(ifsh)*al_key;
      for (iyr=1;iyr<=nyrs_fsh_comp(ifsh);iyr++){
       comp_like_fsh(ifsh) -= nsmpl_fsh(ifsh,iyr)*(oc_fsh(ifsh,iyr)(1,nl_bins) + 0.0001) * 
                              log(ec_fsh(ifsh,iyr)  + 0.0001 ) ;
       }     
     }
     comp_like_fsh(ifsh) -= offset_fsh(ifsh);
   }
   
  if (DebugOut > 0 ) {
    cout << "End Comp_Like" << endl;
    }

  // ============================
FUNCTION Srv_Like
  // ============================

  // Fit to indices (Normal) -------------------------------------------
  if (DebugOut > 0 ) cout << "Begin Srv_Like " << endl;
  surv_like.initialize();
  //surv_net_like.initialize();
  for (int isrv=1;isrv<=nsrv;isrv++)
  for (int iarea=1;iarea<=nareas;iarea++)
    for (int iyr=1;iyr<=nyrs_srv(isrv,iarea);iyr++)
      //surv_like(isrv,iarea) += square(obs_srv(isrv,iarea,iyr) - pred_srv(isrv,iarea,yrs_srv(isrv,iarea,iyr))) / 
      //                             (2.*square(obs_se_srv(isrv,iarea,iyr)));
      surv_like(isrv,iarea) += square(log(obs_srv(isrv,iarea,iyr)) - log(pred_srv(isrv,iarea,yrs_srv(isrv,iarea,iyr)))+0.01) / 
                                   (2.*log(1.0+square(obs_se_srv(isrv,iarea,iyr)/obs_srv(isrv,iarea,iyr))));
  if (DebugOut > 0 ) {
    cout << "pred_srv = " << endl << pred_srv<< endl;
    cout << "obs_srv = " << endl << obs_srv << endl;
    cout << "End Srv_Like " << endl;
    }

  // ============================
FUNCTION Selcofs_Pen 
  // ============================

  selcofs_pen_srv.initialize();
  if (DebugOut > 0 ) cout << "Begin Selcofs_Pen " << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
  {
    if (active(log_selcoffs_srv(isrv)))
    {
      for (int i=1;i<=n_sel_ch_srv(isrv);i++)
      {
        int iyr = yrs_sel_ch_srv(isrv,i) ;
        selcofs_pen_srv(isrv,1) += curv_pen_srv(isrv)*norm2(first_difference(
                                                 first_difference(log_sel_srv(isrv,iyr))));
        if (i>1)
        {
          dvariable var_tmp = square(sel_change_in_srv(isrv,iyr ));
          selcofs_pen_srv(isrv,2)    += .5*norm2( log_sel_srv(isrv,iyr-1) - log_sel_srv(isrv,iyr) ) 
                                   / var_tmp ;
        }
        int nagestmp = nselages_srv(isrv,1);
        for (int j=seldecage;j<=nagestmp;j++)
        {
          dvariable difftmp = log_sel_srv(isrv,iyr,j-1)-log_sel_srv(isrv,iyr,j) ;
          if (difftmp > 0.)
            selcofs_pen_srv(isrv,3)    += .5*square( difftmp ) / seldec_pen_srv(isrv);
        }
        obj_fun            += 20. * square(avgsel_srv(isrv,i));  // To normalize selectivities
      }
    }
  }

  
  // ============================
FUNCTION Fmort_Pen
  // ============================
  dvariable totalN, totalF, TotalG, Temp;
  
  if (DebugOut == 1) cout << "Begin Fmort_Pen" << endl;
  if (active(log_avg_fmort))
   {  
    fpen.initialize();
    for(iarea=1; iarea<=nareas; iarea++)
     {      
      ifsh = iarea;
      if (current_phase()<3) // penalize High F's for beginning phases
       fpen(iarea,1) += 10. * sum(square(Fmort(iarea) - .2));
      else 
       fpen(iarea,1) +=  .001*norm2( Fmort(iarea) - .2);
     }
    for (ifsh = 1; ifsh <= nfsh; ifsh++)
     {
      iarea = areas_fsh(ifsh);
      fpen(iarea,2) += 20*square(mean(fmort_dev(ifsh)) );    
     } 
   }

  // ============================
FUNCTION Compute_priors
  // ============================

  post_priors.initialize();
  //post_priors_srvq.initialize();
  for (isrv=1;isrv<=nsrv;isrv++)
   if (active(log_q_srv(isrv)))
    post_priors_srvq(isrv) += square(q_srv(isrv)-qprior(isrv))/(2*cvqprior(isrv)*cvqprior(isrv)); 
  if (active(steepness))
     post_priors(2) += square(steepness-steepnessprior)/(2*cvsteepnessprior*cvsteepnessprior); 
  for(iarea=1; iarea<=1; iarea++)
   {
    if (active(log_M))
     for (iyr=styr_rec;iyr<=endyr;iyr++)
     post_priors(1) += square(M(iarea,iyr)-mfexp(natmortprior))/(2*cvnatmortprior*cvnatmortprior); 
    if (active(log_sigmar))
     post_priors(3) += square(sigmar(iarea)-sigmarprior(iarea))/(2*cvsigmarprior(iarea)*cvsigmarprior(iarea)); 
   } 

  // ============================
FUNCTION write_mceval
  // ============================

  for(iarea=1;iarea<=nareas;iarea++){
    mceval << "Sp_Biom_mcmc[" << iters << "," << iarea << ",] <- c(";
    for(iyr=styr_sp;iyr<=endyr-1;iyr++)
      mceval << Sp_Biom(iarea,iyr) << ", ";
      mceval << Sp_Biom(iarea,endyr) << ")" << endl;
    }
  for(iarea=1;iarea<=nareas;iarea++){
    mceval << "mod_rec_mcmc[" << iters << "," << iarea << ",] <- c(";
    for(iyr=styr_rec;iyr<=endyr-1;iyr++)
      mceval << mod_rec(iarea,iyr) << ", ";
      mceval << mod_rec(iarea,endyr) << ")" << endl;
    }
  for(iarea=1;iarea<=nareas;iarea++){
    mceval << "M_mcmc[" << iters << "," << iarea << ",] <- c(";
    for(iyr=styr_rec;iyr<=endyr-1;iyr++)
      mceval << M(iarea,iyr)  << ", ";
      mceval << M(iarea,endyr) << ")" << endl;
    }
  for(iarea=1;iarea<=nareas;iarea++){
    mceval << "Bzero_mcmc[" << iters << "," << iarea << "] <- c(";
      mceval << Bzero(iarea) << ")" << endl;
    }
  for(iarea=1;iarea<=nareas;iarea++){
    mceval << "Btot_mcmc[" << iters << "," << iarea << "] <- c(";
      mceval << Btot(iarea) << ")" << endl;
    }
  for(iarea=1;iarea<=nareas;iarea++){
    mceval << "log_Rzero_mcmc[" << iters << "," << iarea << "] <- c(";
      mceval << log_Rzero(iarea) << ")" << endl;
    }
  for(iarea=1;iarea<=nareas;iarea++){
    mceval << "mean_log_rec_mcmc[" << iters << "," << iarea << "] <- c(";
      mceval << mean_log_rec(iarea) << ")" << endl;
    }
  for(iarea=1;iarea<=nareas;iarea++){
    mceval << "log_avg_fmort_mcmc[" << iters << "," << iarea << "] <- c(";
      mceval << log_avg_fmort(iarea) << ")" << endl;
    }

  // MCMC outputs added 7/13/17  

  mceval << "log_Linf_mcmc[" << iters << "] <- c(";
    mceval << log_Linf << ")" << endl;
  mceval << "log_vB_k_mcmc[" << iters << "] <- c(";
    mceval << log_vB_k << ")" << endl;
  mceval << "log_vB_sig_mcmc[" << iters << "] <- c(";
    mceval << log_vB_sig << ")" << endl;
  mceval << "steepness_mcmc[" << iters << "] <- c(";
    mceval << steepness << ")" << endl;    
  for(iarea=1;iarea<=nareas;iarea++){
    mceval << "log_sigmar_mcmc[" << iters << "," << iarea << "] <- c(";
      mceval << log_sigmar(iarea) << ")" << endl;
    }
  for(isrv=1;isrv<=nsrv;isrv++){
    mceval << "log_q_srv_mcmc[" << iters << "," << isrv << "] <- c(";
      mceval << log_q_srv(isrv) << ")" << endl;
    }
  for(isrv=1;isrv<=nsrv;isrv++){
    mceval << "logsel_slope_srv_mcmc[" << iters << "," << isrv << "] <- c(";
      mceval << logsel_slope_srv(isrv) << ")" << endl;
    }
  for(isrv=1;isrv<=nsrv;isrv++){
    mceval << "sel50_srv_mcmc[" << iters << "," << isrv << "] <- c(";
      mceval << sel50_srv(isrv) << ")" << endl;
    }
  for(ifsh=1;ifsh<=nfsh;ifsh++){
    mceval << "log_q_fsh_mcmc[" << iters << "," << ifsh << "] <- c(";
      mceval << log_q_fsh(ifsh) << ")" << endl;
    }
  for(ifsh=1;ifsh<=nfsh;ifsh++){
    mceval << "logsel_slope_fsh_mcmc[" << iters << "," << ifsh << "] <- c(";
      mceval << logsel_slope_fsh(ifsh) << ")" << endl;
    }
  for(ifsh=1;ifsh<=nfsh;ifsh++){
    mceval << "sel50_fsh_mcmc[" << iters << "," << ifsh << "] <- c(";
      mceval << sel50_fsh(ifsh) << ")" << endl;
    }

  for(ifsh=1;ifsh<=nfsh;ifsh++){
  mceval << "fmort_dev_mcmc[" << iters << "," << ifsh << ",] <- c(";
    for(iyr=styr;iyr<=end_datyrs-1;iyr++)
      mceval << fmort_dev(ifsh,iyr) << ", ";
    mceval << fmort_dev(ifsh,end_datyrs) << ")" << endl;
    }


  // ============================
REPORT_SECTION
  // ============================
 
  cout << "starting report" << endl;
  report << " #### INDEX VALUES ##### " << endl;
  report << "nsrv <- " << nsrv << endl;
  report << "nfsh <- " << nfsh << endl;
  report << "styr <- " << styr << endl;
  report << "endyr <- " << endyr << endl;
  report << "end_datyrs <- " << end_datyrs << endl;
  report << "styr_rec <- " << styr_rec << endl;
  report << "styr_sp <- " << styr_sp << endl;
  report << "rec_age <- " << rec_age << endl;
  report << "oldest_age <- " << oldest_age << endl;
  report << "nages <- " << nages << endl;
  report << "nareas <- " << nareas << endl;

  report << endl << "al_key<-list()" << endl;
  report << "al_key <- c(" << endl << al_key(1,1);
  for (age=1;age<=nages;age++){
     if (age == 1)
       istart = 2;
      else
       istart = 1;
      for (int icmp=istart;icmp<=nl_bins;icmp++)
       {
         report << ", " << al_key(age,icmp);        
       }
      }
  report  << endl << ")" << endl;               // al_key
  report << "dim(al_key)<-c(" << nl_bins;
  report << ", " << nages << ")" << endl;
  report << "al_key <- t(al_key)" << endl;
  report << "lbins <- c(" << lbins(1);         // lbins
  for (int icmp=2;icmp<=nl_bins;icmp++)
   {
     report << ", " << lbins(icmp);        
   }
  report << ")" << endl;
  report << "natage <- list()" << endl;         // natage
  for (iarea=1;iarea<=nareas;iarea++)
  {
    report << "natage[[" << iarea << "]] <- list()" << endl;
    report << "natage[[" << iarea << "]] <- c(" << endl << natage(iarea,styr,1);
    for (iyr=styr;iyr<=endyr;iyr++)
      {
       if (iyr == styr)
         istart = 2;
       else
         istart = 1;
       for (age=istart;age<=nages;age++)
         {
          report <<  ", " << natage(iarea,iyr,age);
         }
       }
     report  << endl << ")" << endl;      
     report << "dim(natage[[" << iarea << "]])<-c(" << nages << ", " << nyrs << ")" << endl;
     report << "natage[[" << iarea << "]] <-t(natage[["<< iarea << "]])" << endl;
  }

  report << "n_areasrv <- list()" << endl;         // n_areasrv
  for (isrv=1;isrv<=nsrv;isrv++)
  {
  report << "n_areasrv[[" << isrv << "]] <- list()" << endl;
  for (iarea=1;iarea<=nareas;iarea++)
  {
    report << "n_areasrv[[" << isrv << "]][[" << iarea << "]] <- list()" << endl;
    report << "n_areasrv[[" << isrv << "]][[" << iarea << "]] <- c(" << endl << n_areasrv(isrv,iarea,styr,1);
    for (iyr=styr;iyr<=endyr;iyr++)
      {
       if (iyr == styr)
         istart = 2;
       else
         istart = 1;
       for (age=istart;age<=nages;age++)
         {
          report <<  ", " << n_areasrv(isrv,iarea,iyr,age);
         }
       }
     report  << endl << ")" << endl;      
     report << "dim(n_areasrv[[" << isrv << "]][[" << iarea << "]])<-c(" << nages << ", " << nyrs << ")" << endl;
     report << "n_areasrv[[" << isrv << "]][[" << iarea << "]] <-t(n_areasrv[[" << isrv << "]][["<< iarea << "]])" << endl;
  }
  }


  report << endl << "ec_srv<-list()" << endl;

  for (isrv=1;isrv<=nsrv;isrv++)
  {
  report << "ec_srv[[" << isrv << "]] <- list()" << endl;
  for (iarea=1;iarea<=nareas;iarea++)
    {
     if (nyrs_srv_comp(isrv,iarea) >0)
     {
      report << "ec_srv[[" << isrv << "]][[" << iarea << "]] <- c(" << endl << ec_srv(isrv,iarea,1,1);
      for (iyr=1;iyr<=nyrs_srv_comp(isrv,iarea);iyr++)
        {
          if (iyr == 1)
            istart = 2;
          else
            istart = 1;

          for (int icmp=istart;icmp<=nl_bins;icmp++)
            {
              report << ", " << ec_srv(isrv,iarea,iyr,icmp);        
            }
         }
          report  << endl << ")" << endl;
      }
     }
  }

  for (isrv=1;isrv<=nsrv;isrv++)
   for (iarea=1;iarea<=nareas;iarea++)
    {
     if (nyrs_srv_comp(isrv,iarea) >0)
     {
      report << "dim(ec_srv[[" << isrv << "]][[" << iarea << "]])<-c(" <<   nl_bins;
      report << ", " << nyrs_srv_comp(isrv,iarea) << ")" << endl;
      report << "ec_srv[[" << isrv << "]][[" << iarea << "]]<-t(ec_srv[[" << isrv << "]][[" << iarea << "]])" << endl;
     }
    }
  
                                                    // eac_srv
  report << endl << "eac_srv<-list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
  {
    report << "eac_srv[[" << isrv << "]] <- list()" << endl;
 

    for (iarea=1;iarea<=nareas;iarea++)
    {
     if (nyrs_srv_comp(isrv,iarea) >0)
     {
        report << "eac_srv[[" << isrv << "]][[" << iarea << "]] <- c(" << endl << eac_srv(isrv,iarea,1,1);
        for (iyr=1;iyr<=nyrs_srv_comp(isrv,iarea);iyr++)
          {
            if (iyr == 1)
              istart = 2;
            else
              istart = 1;

            for (int icmp=istart;icmp<=nages;icmp++)
              {
                report << ", " << eac_srv(isrv,iarea,iyr,icmp);        
              }
          }
          report  << endl << ")" << endl;
      }
    }
  }

  for (isrv=1;isrv<=nsrv;isrv++)
   for (iarea=1;iarea<=nareas;iarea++)
    {
     if (nyrs_srv_comp(isrv,iarea) >0)
     {
      report << "dim(eac_srv[[" << isrv << "]][[" << iarea << "]])<-c(" << nages;
      report << ", " << nyrs_srv_comp(isrv,iarea) << ")" << endl;
      report << "eac_srv[[" << isrv << "]][[" << iarea << "]]<-t(eac_srv[[" << isrv << "]][[" << iarea << "]])" << endl;
     }
    }
   
                                                    // oc_srv
  report << endl << "oc_srv<-list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
  {
  report << "oc_srv[[" << isrv << "]] <- list()" << endl;
  for (iarea=1;iarea<=nareas;iarea++)
    {
     if (nyrs_srv_comp(isrv,iarea) >0)
      {
      report << "oc_srv[[" << isrv << "]][[" << iarea << "]] <- c(" << endl <<oc_srv(isrv,iarea,1,1);
      for (iyr=1;iyr<=nyrs_srv_comp(isrv,iarea);iyr++)
        {
          if (iyr == 1)
            istart = 2;
          else
            istart = 1;

          for (int icmp=istart;icmp<=nl_bins;icmp++)
            {
              report << ", " << oc_srv(isrv,iarea,iyr,icmp);        
            }
        }
          report  << endl << ")" << endl;
       }
    }
  }
  for (isrv=1;isrv<=nsrv;isrv++)
  for (iarea=1;iarea<=nareas;iarea++)
    {
     if (nyrs_srv_comp(isrv,iarea) >0)
      {
      report << "dim(oc_srv[[" << isrv << "]][[" << iarea << "]])<-c(" <<   nl_bins;
      report << ", " << nyrs_srv_comp(isrv,iarea) << ")" << endl;
      report << "oc_srv[[" << isrv << "]][[" << iarea << "]]<-t(oc_srv[[" << isrv << "]][[" << iarea << "]])" << endl;
      }
    }
  report << endl;

                                                    // nyrs_srv
  report << endl << "nyrs_srv<-list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
  {
   report << "nyrs_srv[[" << isrv << "]] <- list()" << endl;
   report << "nyrs_srv[[" << isrv << "]][[" << 1 << "]] <- c(" << nyrs_srv(isrv,1);
   for (iarea=2;iarea<=nareas;iarea++)
    {
     report << ", " << nyrs_srv(isrv,iarea);        
    }
    report << ")" << endl;
   }
  report << endl;

                                                    // nyrs_srv_comp
  report << endl << "nyrs_srv_comp<-list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
  {
   report << "nyrs_srv_comp[[" << isrv << "]] <- list()" << endl;
   report << "nyrs_srv_comp[[" << isrv << "]][[" << 1 << "]] <- c(" << nyrs_srv_comp(isrv,1);
   for (iarea=2;iarea<=nareas;iarea++)
    {
      report << ", " << nyrs_srv_comp(isrv,iarea);        
    }
    report << ")" << endl;
   }
  report << endl;
                                                    // yrs_srv_comp
  report << endl << "yrs_srv_comp<-list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
  {
  report << "yrs_srv_comp[[" << isrv << "]] <- list()" << endl;
  for (iarea=1;iarea<=nareas;iarea++)
    {
     if (nyrs_srv_comp(isrv,iarea) >0)
     {
      report << "yrs_srv_comp[[" << isrv << "]][[" << iarea << "]] <- c(" << yrs_srv_comp(isrv,iarea,1);
      for (iyr=2;iyr<=nyrs_srv_comp(isrv,iarea);iyr++)
        {
         report << ", " << yrs_srv_comp(isrv,iarea,iyr);        
        }
          report << ")" << endl;
     }
    }
  report << endl;
  }

  report << endl << "nsmpl_srv<-list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
  {
  report << "nsmpl_srv[[" << isrv << "]] <- list()" << endl;
  for (iarea=1;iarea<=nareas;iarea++)
    {
     if (nyrs_srv_comp(isrv,iarea) >0)
     {
      report << "nsmpl_srv[[" << isrv << "]][[" << iarea << "]] <- c(" << nsmpl_srv(isrv,iarea,1);
      for (iyr=2;iyr<=nyrs_srv_comp(isrv,iarea);iyr++)
        {
         report << ", " << nsmpl_srv(isrv,iarea,iyr);        
        }
          report << ")" << endl;
     }
    }
  }
  report << endl;
                                                    // pred_srv
  report << endl << "pred_srv<-list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
  {
    report << "pred_srv[[" << isrv << "]]<- list()" << endl;
    for (iarea=1;iarea<=nareas;iarea++)
         { 
           report << "pred_srv[[" << isrv << "]][[" << iarea << "]]<- c(" << endl << pred_srv(isrv,iarea,styr);
           for (iyr=styr+1;iyr<=endyr;iyr++)
             {
               report << ", " << pred_srv(isrv,iarea,iyr);
             }
           report << ")" << endl;
          }        
  }
                                                   // yrs_srv
  // init_3darray  yrs_srv(1,nsrv,1,nareas,1,nyrs_srv)
  report << endl << "yrs_srv <- list()" << endl;
   for (isrv=1;isrv<=nsrv;isrv++)
   {
     report << "yrs_srv[[" << isrv << "]] <- list()" << endl;
     for (iarea=1;iarea<=nareas;iarea++)
      {
       if (nyrs_srv(isrv,iarea) >0)
       {
        report << "yrs_srv[[" << isrv << "]][[" << iarea << "]] <- c(" << yrs_srv(isrv,iarea,1);
        for(iyr=2;iyr<=nyrs_srv(isrv,iarea);iyr++)
         {
          report << ", " << yrs_srv(isrv,iarea,iyr);
         }
         report   << ")" << endl;
        }
      }
   }
  report << endl;
                                                    // obs_srv
  report << "obs_srv <- list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
  {
   report << "obs_srv[[" << isrv << "]] <- list()" << endl;
   for (iarea=1;iarea<=nareas;iarea++)
    {
     if (nyrs_srv(isrv,iarea) >0)
     {
      report << "obs_srv[[" << isrv << "]][[" << iarea << "]] <- c(" << obs_srv(isrv,iarea,1);
      for(iyr=2;iyr<=nyrs_srv(isrv,iarea);iyr++)
       {
        report << ", " << obs_srv(isrv,iarea,iyr);
       }
          report   << ")" << endl;
      }
     }
    }
  report << endl;
                                                   //obs_se_srv
  report << "obs_se_srv <- list()" << endl;
   for (isrv=1;isrv<=nsrv;isrv++)
  {
   report << "obs_se_srv[[" << isrv << "]] <- list()" << endl;
   for (iarea=1;iarea<=nareas;iarea++)
    {
     if (nyrs_srv(isrv,iarea) >0)
     {
      report << "obs_se_srv[[" << isrv << "]][[" << iarea << "]] <- c(" << obs_se_srv(isrv,iarea,1);
      for(iyr=2;iyr<=nyrs_srv(isrv,iarea);iyr++)
       {
        report << ", " << obs_se_srv(isrv,iarea,iyr);
       }
       report  << ")" << endl;
      }
    }
  }
  report << endl;


  // FISHERIES
  cout << "starting fisheries part of report" << endl;
  report << "ec_fsh <- list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
    {
     if (nyrs_fsh_comp(ifsh) >0)
     {
      report << "ec_fsh[[" << ifsh << "]] <- c(" << endl << ec_fsh(ifsh,1,1);
      for (iyr=1;iyr<=nyrs_fsh_comp(ifsh);iyr++)
        {
          if (iyr == 1)
            istart = 2;
          else
            istart = 1;

          for (int icmp=istart;icmp<=nl_bins;icmp++)
            {
              report << ", " << ec_fsh(ifsh,iyr,icmp);        
            }
         }
          report  << endl << ")" << endl;
      }
     }
  }

  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
     if (nyrs_fsh_comp(ifsh) >0)
     {
      report << "dim(ec_fsh[[" << ifsh << "]])<-c(" <<   nl_bins;
      report << ", " << nyrs_fsh_comp(ifsh) << ")" << endl;
      report << "ec_fsh[[" << ifsh << "]]<-t(ec_fsh[[" << ifsh << "]])" << endl;
     }
    }
  
                                                    // eac_fsh
  report << endl << "eac_fsh<-list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
    //report << "eac_fsh[[" << ifsh << "]] <- list()" << endl;
     {
     if (nyrs_fsh_comp(ifsh) >0)
     {
        report << "eac_fsh[[" << ifsh << "]] <- c(" << endl << eac_fsh(ifsh,1,1);
        for (iyr=1;iyr<=nyrs_fsh_comp(ifsh);iyr++)
          {
            if (iyr == 1)
              istart = 2;
            else
              istart = 1;

            for (int icmp=istart;icmp<=nages;icmp++)
              {
                report << ", " << eac_fsh(ifsh,iyr,icmp);        
              }
          }
          report  << endl << ")" << endl;
      }
    }
  }
  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
     if (nyrs_fsh_comp(ifsh) >0)
     {
      report << "dim(eac_fsh[[" << ifsh << "]])<-c(" << nages;
      report << ", " << nyrs_fsh_comp(ifsh) << ")" << endl;
      report << "eac_fsh[[" << ifsh << "]]<-t(eac_fsh[[" << ifsh << "]])" << endl;
     }
    }
                                                    // oc_fsh
  report << endl << "oc_fsh<-list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
 // report << "oc_fsh[[" << ifsh << "]] <- list()" << endl;
    {
     if (nyrs_fsh_comp(ifsh) >0)
      {
      report << "oc_fsh[[" << ifsh << "]] <- c(" << endl <<oc_fsh(ifsh,1,1);
      for (iyr=1;iyr<=nyrs_fsh_comp(ifsh);iyr++)
        {
          if (iyr == 1)
            istart = 2;
          else
            istart = 1;

          for (int icmp=istart;icmp<=nl_bins;icmp++)
            {
              report << ", " << oc_fsh(ifsh,iyr,icmp);        
            }
        }
          report  << endl << ")" << endl;
       }
    }
  }
  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
     if (nyrs_fsh_comp(ifsh) >0)
      {
      report << "dim(oc_fsh[[" << ifsh << "]])<-c(" <<   nl_bins;
      report << ", " << nyrs_fsh_comp(ifsh) << ")" << endl;
      report << "oc_fsh[[" << ifsh << "]]<-t(oc_fsh[[" << ifsh << "]])" << endl;
      }
    }
  report << endl;

                                                    // nyrs_fsh
  report << endl << "nyrs_fsh<-vector()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
   //report << "nyrs_fsh[[" << ifsh << "]] <- list()" << endl;
   report << "nyrs_fsh <- c(" << nyrs_fsh(1);
   if (nfsh >1) 
    {
     report << ", " << nyrs_fsh(ifsh);        
    }
   }
  report << ")" << endl;
  report << endl;

                                                    // nyrs_fsh_comp
  report << endl << "nyrs_fsh_comp<-vector()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
   report << "nyrs_fsh_comp <- c(";
    {
      report  << nyrs_fsh_comp(ifsh);
      if(ifsh < nfsh) report << ", ";        
    }
    report << ")" << endl;
   }
  report << endl;
                                                    // yrs_fsh_comp
  report << endl << "yrs_fsh_comp<-list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
  report << "yrs_fsh_comp[[" << ifsh << "]] <- vector()" << endl;
    {
     if (nyrs_fsh_comp(ifsh) >0)
     {
      report << "yrs_fsh_comp[[" << ifsh << "]] <- c(" << yrs_fsh_comp(ifsh,1);
      for (iyr=2;iyr<=nyrs_fsh_comp(ifsh);iyr++)
        {
         report << ", " << yrs_fsh_comp(ifsh,iyr);        
        }
          report << ")" << endl;
     }
    }
  report << endl;
  }

  report << endl << "nsmpl_fsh<-list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
  report << "nsmpl_fsh[[" << ifsh << "]] <- vector()" << endl;
    {
     if (nyrs_fsh_comp(ifsh) >0)
     {
      report << "nsmpl_fsh[[" << ifsh << "]] <- c(" << nsmpl_fsh(ifsh,1);
      for (iyr=2;iyr<=nyrs_fsh_comp(ifsh);iyr++)
        {
         report << ", " << nsmpl_fsh(ifsh,iyr);        
        }
          report << ")" << endl;
     }
    }
  }
  report << endl;

  // init_3darray  yrs_fsh(1,nfsh,1,nareas,1,nyrs_fsh)
  report << endl << "yrs_fsh <- list()" << endl;
   for (ifsh=1;ifsh<=nfsh;ifsh++)
   {
     report << "yrs_fsh[[" << ifsh << "]] <- vector()" << endl;
      {
       if (nyrs_fsh(ifsh) >0)
       {
        report << "yrs_fsh[[" << ifsh << "]] <- c(" << yrs_fsh(ifsh,1);
        for(iyr=2;iyr<=nyrs_fsh(ifsh);iyr++)
         {
          report << ", " << yrs_fsh(ifsh,iyr);
         }
         report   << ")" << endl;
        }
      }
   }
  report << endl;
                                                    // obs_catch
  report << "obs_catch <- list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
   report << "obs_catch[[" << ifsh << "]] <- vector()" << endl;
    {
     //if (nyrs_fsh(ifsh) >0)
     {
      report << "obs_catch[[" << ifsh << "]] <- c(" << obs_catch(ifsh,styr);
      for(iyr=styr+1;iyr<=endyr;iyr++)
       {
        report << ", " << obs_catch(ifsh,iyr);
       }
          report   << ")" << endl;  
      }
     }
    }
  report << endl;
                                                   //obs_se_catch
  report << "obs_se_catch <- vector()" << endl;
      report << "obs_se_catch[1] <- c(" << obs_se_catch(1);
      if(nareas > 1)
      for(iarea =2;iarea<=nareas;iarea++)
       {
        report << ", " << obs_se_catch(iarea);
       }
       report  << ")" << endl;
  report << endl;

  report << "c.gamma <- " << gamma << endl;

  cout << "fisheries predictions report" << endl;
                                                     // pred_catch
  //3darray            pred_catch(1,lnfsh,1,lnareas,styr,endyr)
  report << "pred_catch <- list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
   report << "pred_catch[[" << ifsh << "]] <- vector()" << endl;
    {
     //if (nyrs_fsh(ifsh) >0)
     {
      report << "pred_catch[[" << ifsh << "]] <- c(" << pred_catch(ifsh,styr);
      for(iyr=styr+1;iyr<=endyr;iyr++)
       {
        report << ", " << pred_catch(ifsh,iyr);
       }
          report   << ")" << endl;  
      }
     }
    }
  report << endl;

  //init_matrix        log_avg_fmort(1,nareas,1,nfsh,phase_fmort)
  report << "log_avg_fmort <- vector()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
   //report << "log_avg_fmort[[" << ifsh << "]] <- list()" << endl;
    {
      report << "log_avg_fmort <- c(" << log_avg_fmort(ifsh);
      if(nfsh > 1)
        for(ifsh=2;ifsh<=nfsh;ifsh++)
          report << ", " << log_avg_fmort(ifsh);
      report << ")" << endl;
     }
    }
  report << endl;
  
  report << "fmort_dev <- list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
   //report << "fmort_dev[[" << ifsh << "]] <- vector()" << endl;
    {
     //if (nyrs_fsh(ifsh) >0)
     {
      report << "fmort_dev[[" << ifsh << "]] <- c(" << fmort_dev(ifsh,styr);
      for(iyr=styr+1;iyr<=end_datyrs;iyr++)
       {
        report << ", " << fmort_dev(ifsh,iyr);
       }
          report   << ")" << endl;  
      }
     }
    }
  report << endl;

  report << "Fmort <- list()" << endl;
  for (iarea=1;iarea<=nareas;iarea++)
    {
     //if (nyrs_fsh(iarea) >0)
     {
      report << "Fmort[[" << iarea << "]] <- c(" << Fmort(iarea,styr);
      for(iyr=styr+1;iyr<=end_datyrs;iyr++)
       {
        report << ", " << Fmort(iarea,iyr);
       }
          report   << ")" << endl;  
      }
     }
  report << endl;

  if (catch_opt == 4)
   {
    //   dvar_matrix fmort_calc(1,nfsh,styr,end_datyrs);// F  for catch_opt == 4
    report << "fmort_calc <- list()" << endl;
    for (ifsh=1;ifsh<=nareas;ifsh++)
    {
      report << "fmort_calc[[" << ifsh << "]] <- c(" << fmort_calc(ifsh,styr);
      for(iyr=styr+1;iyr<=end_datyrs;iyr++)
       {
        report << ", " << fmort_calc(ifsh,iyr);
       }
       report   << ")" << endl;  
      }
    report << endl;
    }

  report << "catage <- list()" << endl;         // catage
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
    //report << "catage[[" << ifsh << "]] <- list()" << endl;
    report << "catage[[" << ifsh <<"]] <- c(" << endl << catage(ifsh,styr,1);
    for (iyr=styr;iyr<=endyr;iyr++)
      {
       if (iyr == styr)
         istart = 2;
       else
         istart = 1;
       for (age=istart;age<=nages;age++)
         {
          report <<  ", " << catage(ifsh,iyr,age);
         }
       }
     report  << endl << ")" << endl;      
     report << "dim(catage[[" << ifsh << "]])<-c(" << nages << ", " << nyrs << ")" << endl;
     report << "catage[[" << ifsh <<"]] <-t(catage[[" << ifsh <<"]])" << endl;
  }
  report << endl;

  report << "catage_fut <- list()" << endl;         // catage_fut
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
    //report << "catage_fut[[" << ifsh << "]] <- list()" << endl;
    report << "catage_fut[[" << ifsh <<"]] <- c(" << endl << catage_fut(ifsh,end_datyrs+1,1);
    for (iyr=end_datyrs+1;iyr<=endyr;iyr++)
      {
       if (iyr == end_datyrs+1)
         istart = 2;
       else
         istart = 1;
       for (age=istart;age<=nages;age++)
         {
          report <<  ", " << catage_fut(ifsh,iyr,age);
         }
       }
     report  << endl << ")" << endl;      
     report << "dim(catage_fut[[" << ifsh << "]])<-c(" << nages << ", " << nproj_yrs << ")" << endl;
     report << "catage_fut[[" << ifsh <<"]] <-t(catage_fut[[" << ifsh <<"]])" << endl;
  }
  report << endl;

  report << endl << "effort<-list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
    report << "effort[[" << ifsh << "]] <- vector()" << endl;
      report << "effort[[" << ifsh << "]] <- c("  << effort(ifsh,styr);
      for (iyr=styr+1;iyr<=endyr;iyr++)
        {
         report << ", " << effort(ifsh,iyr);        
        }
      report  << endl << ")" << endl;
    }
 
  report << "F <- list()" << endl;         // F
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
    report << "F[[" << ifsh <<"]] <- c(" << endl << F(ifsh,styr,1);
    for (iyr=styr;iyr<=end_datyrs;iyr++)
      {
       if (iyr == styr)
         istart = 2;
       else
         istart = 1;
       for (age=istart;age<=nages;age++)
         {
          report <<  ", " << F(ifsh,iyr,age);
         }
       }
     report  << endl << ")" << endl;      
     report << "dim(F[[" << ifsh << "]])<-c(" << nages << ", " << ndatyrs << ")" << endl;
     report << "F[[" << ifsh <<"]] <-t(F[[" << ifsh <<"]])" << endl;
  }

                                                     // Sp_Biom
  report << endl << "Sp_Biom<-list()" << endl;
    for(iarea=1;iarea<=nareas;iarea++)
    {
      report << "Sp_Biom[[" << iarea << "]] <- c(" << endl << Sp_Biom(iarea,styr);
      for (iyr=styr_sp+1; iyr<=endyr; iyr++)
       {
        report << ", " << Sp_Biom(iarea,iyr);
       }
       report << ")" << endl;
    }

                                                    // pred_rec
  report << endl << "pred_rec<-list()" << endl;
    for(iarea=1;iarea<=nareas;iarea++)
    {
      iyr = styr_rec;
      report << "pred_rec[[" << iarea << "]] <- c(" << endl << pred_rec(iarea,iyr);
      for (iyr=styr_rec+1; iyr<=endyr; iyr++)
       {
         report << ", " << pred_rec(iarea,iyr);
       }
       report << ")" << endl;
    }

                                                    // mod_rec
  report << endl << "mod_rec<-list()" << endl;
    for(iarea=1;iarea<=nareas;iarea++)
    {
      iyr = styr_rec;
      report << "mod_rec[[" << iarea << "]] <- c(" << endl << mod_rec(iarea,iyr);
      for (iyr=styr_rec+1; iyr<=endyr; iyr++)
       {
         report << ", " << mod_rec(iarea,iyr);
       }
       report << ")" << endl;
    }
                                                    // sel_srv
  report << endl << "sel_srv<-list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
         { 
           msrv = msrv_id(isrv);
           report << "sel_srv[[" << isrv << "]]<- c(" << endl << sel_srv(msrv,styr,1);
           for (iyr=styr;iyr<=endyr;iyr++)
             {
               if (iyr == styr)
                 istart = 2;
               else
                 istart = 1;

               for (int icmp=istart;icmp<=nages;icmp++)
                 {
                   report << ", " << sel_srv(msrv,iyr,icmp);        
                 }
              }
          report  << endl << ")" << endl;
         }

  for (isrv=1;isrv<=nsrv;isrv++)
    {
      report << "dim(sel_srv[[" << isrv << "]])<-c(" << nages;
      report << ", " << nyrs << ")" << endl;
      report << "sel_srv[[" << isrv << "]]<-t(sel_srv[[" << isrv << "]])" << endl;
    }

                                                    // log_sel_fsh
  report << endl << "log_sel_fsh<-list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
         { 
           report << "log_sel_fsh[[" << ifsh << "]]<- c(" << endl << log_sel_fsh(ifsh,styr,1);
           for (iyr=styr;iyr<=endyr;iyr++)
             {
               if (iyr == styr)
                 istart = 2;
               else
                 istart = 1;

               for (int icmp=istart;icmp<=nages;icmp++)
                 {
                   report << ", " << log_sel_fsh(ifsh,iyr,icmp);        
                 }
              }
          report  << endl << ")" << endl;
         }

  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      report << "dim(log_sel_fsh[[" << ifsh << "]])<-c(" << nages;
      report << ", " << nyrs << ")" << endl;
      report << "log_sel_fsh[[" << ifsh << "]]<-t(log_sel_fsh[[" << ifsh << "]])" << endl;
    }

  cout << "report biological parameters" << endl;
  report << " #### ESTIMATED BIOLOGICAL PARAMETERS #### " << endl;
  //report << "M <- " << M << endl;
  // M(1,nareas,styr_rec,endyr)
  //report << "M<-c(" << M(1);
  //for (iarea=2;iarea<=nareas;iarea++)
   //report << ", " << M(iarea);
  //report << ")" << endl; 
  report << "M <- array(dim=c(" << nareas << ", " << endyr-styr_rec+1 << "))"<< endl;
    for (iarea=1;iarea<=nareas;iarea++)
     {
      report << "M[" << iarea << ",] <- c(" ;
      for (iyr=styr_rec;iyr<=endyr;iyr++)
      {
       report << M(iarea,iyr);
       if (iyr < endyr)  report << ", ";
       else if(iyr==endyr) report << ")" << endl; 
      } 
     }

  report << "mean_log_rec <- c(" << mean_log_rec(1);
  for (iarea=2;iarea<=nareas;iarea++)
    report << ", " << mean_log_rec(iarea);
  report << ")" << endl;

  report << "steepness<- c(" << steepness << ")" << endl;

  report << "log_Rzero<- c(" << log_Rzero(1);
  for (iarea=2;iarea<=nareas;iarea++)
     report << ", " << log_Rzero(iarea);
  report << ")" << endl;

  report << "Bzero<- c(" << Bzero(1);
  for (iarea=2;iarea<=nareas;iarea++)
     report << ", " << Bzero(iarea);
  report << ")" << endl; 

  report << "Btot<- c(" << Btot(1);
  for (iarea=2;iarea<=nareas;iarea++)
     report << ", " << Btot(iarea);
  report << ")" << endl; 


  report << "rec_dev <- array(dim=c(" << nareas << ", " << end_datyrs-styr_rec+1 << "))"<< endl;
    for (iarea=1;iarea<=nareas;iarea++)
     {
      report << "rec_dev[" << iarea << ",] <- c(" ;
      for (iyr=styr_rec;iyr<=end_datyrs;iyr++)
      {
       report << rec_dev(iarea,iyr);
       if (iyr < end_datyrs)  report << ", ";
       else if(iyr==end_datyrs) report << ")" << endl; 
      } 
     }
  report << "nproj_yrs <- " << nproj_yrs << endl;

  report << "log_sigmar<- c(" << log_sigmar(1);
  for (iarea=2;iarea<=nareas;iarea++)
    report << ", " << log_sigmar(iarea);
  report << ")" << endl; 
  report << endl; 
                                                    // wt_len
  report << "wt_len <- c(" << wt_len(1);
  for (ilen=2;ilen<=nl_bins;ilen++)
    report << ", " << wt_len(ilen);
  report << ")" << endl;
                                                    // wt_age
  report << "wt_age <- c(" << wt_age(1);
  for (age=2;age<=nages;age++)
    report << ", " << wt_age(age);
  report << ")" << endl << endl;
 
  report << " ## movement parameters ## " << endl;
  if(nareas > 1){                                                      // nu
  report << "nu <- array(dim=c(" << nareas << ", " << nareas-1 << ", " << endyr-styr_rec+1 << "))"<< endl;
    for (iarea=1;iarea<=nareas;iarea++)
     {
      for(iyr=styr_rec;iyr<= endyr;iyr++)
       {
        report << "nu[" << iarea << ",," << endyr-iyr+1 << "] <- c(" ;
        for (jarea=1;jarea<=(nareas-1);jarea++)
        {
         report << nu(iarea+(jarea-1),iyr);
         if (jarea < (nareas-1))  
           report << ", ";
         else 
           report << ")" << endl; 
        } 
       }
      }
                                                     // mean_log_move
  report << "mean_log_move <- array(dim=c(" << nareas << ", " << nareas-1 << "))"<< endl;
    for (iarea=1;iarea<=nareas;iarea++)
       {
        report << "mean_log_move[" << iarea << "," << "] <- c(" ;
        for (jarea=1;jarea<=(nareas-1);jarea++)
        {
         report << mean_log_move(iarea+jarea-1);
         if (jarea < (nareas-1))  report << ", ";
         else if(jarea==(nareas-1)) report << ")" << endl; 
        } 
       }
                                                    // log_move_dev
  report << "log_move_dev <- array(dim=c(" << nareas << ", " << nareas-1 << ", " << endyr-styr_rec+1 << "))"<< endl;
    for (iarea=1;iarea<=nareas;iarea++)
      for(iyr=styr_rec;iyr<= endyr;iyr++)
       {
        report << "log_move_dev[" << iarea << ",," << iyr-styr_rec+1 << "] <- c(" ;
        for (jarea=1;jarea<=(nareas-1);jarea++)
        {
         report << log_move_dev(iarea+jarea-1,iyr);
         if (jarea < (nareas-1))  report << ", ";
         else if(jarea==(nareas-1)) report << ")" << endl; 
        } 
       }
  }
  if(nareas == 1){                                                      // nu
  report << "nu <- vector(length = "  << endyr-styr_rec+1 << ")" << endl;
  report << endl; 
   report << "mean_log_move <- " << mean_log_move << endl;
                                                    // log_move_dev
  report << "log_move_dev <- array(dim=c(" << nareas << ", " << endyr-styr_rec+1 << "))" << endl;
  for(iarea=1;iarea<=nareas;iarea++)
   {
    report << "log_move_dev[" << iarea << ",]<- c(" << log_move_dev(iarea,styr_rec);
      for(iyr=styr_rec+1;iyr<=endyr;iyr++)
       report << ", " << log_move_dev(iarea,iyr);
    report << ")" << endl; 
   }
  report << endl; 
  }
  
  report << "immigration <- list()" << endl;         // immigration
  for (iarea=1;iarea<=nareas;iarea++)
  {
    report << "immigration[[" << iarea << "]] <- list()" << endl;
    report << "immigration[[" << iarea << "]] <- c(" << endl << immigration(iarea,styr,1);
    for (iyr=styr_rec;iyr<=endyr;iyr++)
      {
       if (iyr == styr)
         istart = 2;
       else
         istart = 1;
       for (age=istart;age<=nages;age++)
         {
          report <<  ", " << immigration(iarea,iyr,age);
         }
       }
     report  << endl << ")" << endl;      
     report << "dim(immigration[[" << iarea << "]])<-c(" << nages << ", " << endyr-styr_rec+1 << ")" << endl;
     report << "immigration[[" << iarea << "]] <-t(immigration[["<< iarea << "]])" << endl;
  }

  //report << "emig_rate <- immig_rate <- array(dim=c(" << nareas << ", " << endyr-styr_rec+1 << "))" << endl;
  report << "emig_rate <- array(dim=c(" << nareas << ", " << endyr-styr_rec+1 << "))" << endl;
  for(iarea=1;iarea<=nareas;iarea++)
   {
    report << "emig_rate[" << iarea << ",]<- c(" << emig_rate(iarea,styr_rec);
      for(iyr=styr_rec+1;iyr<=endyr;iyr++)
       report << ", " << emig_rate(iarea,iyr);
    report << ")" << endl; 
   }
  report << endl;

  // init_number_vector log_q_srv(1,nsrv,phase_q_srv)
  report << "q_srv <- c(" << mfexp(log_q_srv[1]);
    for (isrv=2;isrv<=nsrv;isrv++)
      report << ", " <<  mfexp(log_q_srv[isrv]);
  report << ")"<< endl;

  report << "q_fsh <- c(" << mfexp(log_q_fsh[1]);
    for (ifsh=2;ifsh<=nfsh;ifsh++)
      report << ", " << mfexp(log_q_fsh[ifsh]);
  report << ")"<< endl;

  report << "influx<- c(" << influx(1);
  for (iarea=2;iarea<=nareas;iarea++)
    report << ", " << influx(iarea);
  report << ")" << endl; 
  report << endl; 
  
  report << "Linf <- " << Linf << endl;
  report << "vB_k <- " << vB_k << endl;
  report << "vB_sig <- " << vB_sig << endl;

  report << "catch_opt <- " << catch_opt << endl;
  report << "move_pen <- " << move_pen << endl;

  report << "catch_like <- array(dim=c(nfsh,nareas))" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
  {
    report << "catch_like[" << ifsh << ",] <-c(" << catch_like(ifsh,1);
    for (iarea=2;iarea<=nareas;iarea++)
     report << ", " << catch_like(ifsh,iarea);
    report << ")" << endl; 
  }
  report << "move_pen <- " << move_pen << endl;
  report << "surv_like <- array(dim=c(nsrv,nareas))" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
  {
    report << "surv_like[" << isrv << ",] <-c(" << surv_like(isrv,1);
    for (iarea=2;iarea<=nareas;iarea++)
     report << ", " << surv_like(isrv,iarea);
    report << ")" << endl; 
  }
  report << "comp_like_srv <-c(" << comp_like_srv(1);
  for (isrv=2;isrv<=nsrv;isrv++)
   report << ", " << comp_like_srv(isrv);
  report << ")" << endl; 
  report << "comp_like_fsh <-c(" << comp_like_fsh(1);
  for (iarea=2;iarea<=nareas;iarea++)
   report << ", " << comp_like_fsh(iarea);
  report << ")" << endl; 
  report << "selcofs_pen_srv <- " <<  sum(selcofs_pen_srv) << endl;
  report << "rec_pen <-c(" << rec_pen(1);
  for (iarea=2;iarea<=4;iarea++)
   report << ", " << rec_pen(iarea);
  report << ")" << endl; 
  report << "CrashPen <- " << CrashPen << endl;
  //report << "selcofs_pen_fsh <- " <<  sum(selcofs_pen_fsh) << endl;
  report << "fpen <-" << sum(fpen) << endl;
  report << "post_priors_srvq <- " << sum(post_priors_srvq) << endl;
  report << "post_priors <- " << sum(post_priors) << endl;

  report << "alpha<- c(" << alpha(1);
  for (iarea=2;iarea<=nareas;iarea++)
     report << ", " << alpha(iarea);
  report << ")" << endl; 
  report << "beta<- c(" << beta(1);
  for (iarea=2;iarea<=nareas;iarea++)
     report << ", " << beta(iarea);
  report << ")" << endl; 

  report << "obj_comps <- c(" << endl;
  for(i=1;i<=11;i++)
    report << obj_comps(i) << ",";
  report << obj_comps(12) << ")" << endl;

  report << "surv_like_sum <- " << sum(surv_like) << endl;
  report << "comp_like_srv_sum <- " << sum(comp_like_srv) << endl;
  report << "max_g <- " << objective_function_value::pobjfun->gmax << endl;
  report << "obj_val <- " << obj_fun << endl;

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(6263);  // replaced 1000 with 6263
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(20000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(15000000);
  arrmblsize=500000000;

