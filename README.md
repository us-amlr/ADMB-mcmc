# ADMB-mcmc
Example MCMC using ADMB code and krill data from Kinzey et al., 2018

This repository provides AD Model Builder (ADMB) code ('model/krill.tpl') and data ('model/krill.dat') for configuration XVI in Kinzey et al., 2018 ('KinzeyWattersReiss_2018.pdf' in this repository). It also contains an executable ('model/krill.exe') that was precompiled using ADMB version 12.3 (version 11.1 for the original paper) and then compiled with MinGW.

To reproduce the model, download 'krill.dat', 'krill.exe' (or compile from 'krill.tpl'), 'krill.pin' (starting values for parameter estimation) and 'mcmc.bat' (ADMB console commands to save 5000 MCMC samples from 5 million iterations) from the 'models' directory in this repository to a local directory and double-click the batch file. After about three hours the 19 GB, R-readable 'mceval.dat' containing 5000 MCMC samples will be produced (the file size and run time can be reduced by decreasing the numbers in 'mcmc.bat'). Additional output files are the maximum likelihood parameter estimates ('krill.par'), their standard deviations ('krill.std'), and an R-readable report file ('krill.rep').

To obtain only the maximum likelihood estimates, use 'krill -ind krill.dat >> out' from a console window. The 'out' command redirects output from the computer screen to a file and so runs faster. This file can be thrown away when the MCMC sampling is finished.

Figures 3c and 4c in Kinzey et al., 2018 illustrate the MCMC samples for numbers of recruits and spawning biomass, respectively ('mcmc_plot_Figs4c_3c.pdf' in the main directory). These can be reproduced from the 'mceval.dat' file using 'mcmc_plot_Figs4c_3c.r'.

There are two diffences in the 'krill.tpl' in this repository and the 'ADMB-krill-spawner-recruit' repository that is also on the AMLR GitHub page. These are in how the 'RecTmp' variable in the case 3 portion of the SRecruit function, and how the survey likelihood (variable 'surv_like'), are defined. These differences suited the different objectives of the two papers.
