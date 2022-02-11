# ADMB-mcmc
Example MCMC using ADMB code and krill data from Kinzey et al., 2018

This repository provides AD Model Builder (ADMB) code ('model/krill.tpl') and data ('model/krill.dat') for configuration XVI in Kinzey et al., 2018 ('KinzeyWattersReiss_2018.pdf' in this repository). It also contains an executable ('model/krill.exe') that was precompiled using ADMB version 12.3 (version 11.1 for the original paper) and then compiled with MinGW.

To reproduce the model, download 'krill.dat', 'krill.exe' (or compile from 'krill.tpl'), 'krill.pin' (starting values for parameter estimation) and 'mcmc.bat' (ADMB console commands to save 5000 MCMC samples from 5 million iterations) in this repository to a local directory. Either double-click the batch file or enter its commands manually from the console. After about three hours the 19 GB, R-readable 'mceval.dat' containing 5000 MCMC samples will be produced (the file size and run time can be reduced by decreasing the numbers in 'mcmc.bat'). Additional output files are the maximum likelihood parameter estimates ('krill.par'), their standard deviations ('krill.std'), and an R-readable report file ('krill.rep').

To obtain only the maximum likelihood estimates without the 'mceval.dat' file, use 'krill -ind krill.dat >> out' from a console window. The 'out' command redirects output from the computer screen to a file and so runs faster. This file can be thrown away when the run is finished.

Figures 3c and 4c in Kinzey et al., 2018 illustrate the MCMC samples for numbers of recruits and spawning biomass, respectively ('mcmc_plot_Figs4c_3c.pdf' in the main directory). These can be reproduced from the 'mceval.dat' file using 'mcmc_plot_Figs4c_3c.r'.

There are two diffences in the 'krill.tpl' in this repository and the 'ADMB-krill-spawner-recruit' repository that is also on the AMLR GitHub page. These are in how the 'RecTmp' variable in the case 3 portion of the SRecruit function, and how the survey likelihood (variable 'surv_like'), are defined. These differences suited the different objectives of the two papers.

# Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
