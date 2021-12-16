rm(list=ls())
#setwd()
iters <- 5000 # number saved mcmc samples (-mcmc/-mcsave from 'mcmc.bat')
path.in <- getwd()
path.out <- paste(path.in,'/Fig4c_mcmc/',sep='')
config.nm <- c(  
               '/model/'
              )
if(!dir.exists('Fig4c_mcmc')){
  dir.create('Fig4c_mcmc')
  }
est.nms <- c('Sp_Biom','mod_rec'
	    )
plot.nms <- c('log Spawning Biomass (tonnes)', 'log Recruit Abundance (N)', 'M')
for(i.cfg in 1:length(config.nm)){
  source(paste(path.in,config.nm[i.cfg],'krill.rep',sep=''))
  yrs <- styr_rec:endyr
  nyrs.2 <- endyr-styr_rec+2
  Sp_Biom_mcmc <-  array(dim=c(iters,nareas,nyrs.2))
  M_mcmc <- mod_rec_mcmc <- array(dim=c(iters,nareas,(nyrs.2-1)))
    Bzero_mcmc <- Btot_mcmc <- log_Rzero_mcmc <- mean_log_rec_mcmc <- log_avg_fmort_mcmc <- 
                  array(dim=c(iters,nareas))
    log_Linf_mcmc <- log_vB_k_mcmc <- log_vB_sig_mcmc <-
           steepness_mcmc <- vector()
    log_sigmar_mcmc <- array(dim=c(iters,nareas))
    log_q_srv_mcmc <- logsel_slope_srv_mcmc <- sel50_srv_mcmc <-
                      array(dim=c(iters,nsrv))
    log_q_fsh_mcmc <- logsel_slope_fsh_mcmc <- sel50_fsh_mcmc <-
                      array(dim=c(iters,nfsh))
    fmort_dev_mcmc <- array(dim=c(iters,nareas,length(styr:end_datyrs))) 

  source(paste(path.in,config.nm[i.cfg],'mceval.dat',sep=''))

  for(i.area in 1:nareas){
    plt.name <- paste(path.out,'PLOTS_mcmc_SpB_rec_',i.cfg,'.pdf',sep='')
    pdf(file = plt.name,width=14,height=5)
    par(mfrow=c(1,3),cex=1,oma=c(2,2,2,0))
  for(i.est in 1:length(est.nms)){
    est_mcmc <- get(paste(est.nms[i.est],'_mcmc',sep=''))
    est_mcmc <- log(est_mcmc[,i.area,])
    x.yrs <- styr_rec:end_datyrs
   { 
    ifelse(
      i.est == 1,
      est_mcmc <- est_mcmc[,2:47], # styr_rec:end_datyrs
      est_mcmc <- est_mcmc[,1:46]
      )
    y.lim <- c(min(apply(est_mcmc,2,quantile,na.rm=T)),
           max(apply(est_mcmc,2,quantile,na.rm=T)))
    plot(x.yrs,apply(est_mcmc,2,quantile,na.rm=T)[1,],
      type='l',ylim=y.lim,ylab='',
      xlab = '',main = plot.nms[i.est])
    polygon(c(x.yrs, rev(x.yrs)), 
      c(apply(est_mcmc,2,quantile,na.rm=T)[2,], 
      rev(apply(est_mcmc,2,quantile,na.rm=T)[1,])),
      col = 'grey92', border = NA)
    polygon(c(x.yrs, rev(x.yrs)), 
      c(apply(est_mcmc,2,quantile,na.rm=T)[3,], 
      rev(apply(est_mcmc,2,quantile,na.rm=T)[2,])),
      col = 'light grey', border = NA)
    polygon(c(x.yrs, rev(x.yrs)), 
      c(apply(est_mcmc,2,quantile,na.rm=T)[5,], 
      rev(apply(est_mcmc,2,quantile,na.rm=T)[4,])),
      col = 'grey92', border = NA)
    polygon(c(x.yrs, rev(x.yrs)), 
      c(apply(est_mcmc,2,quantile,na.rm=T)[4,], 
      rev(apply(est_mcmc,2,quantile,na.rm=T)[3,])),
      col = 'light grey', border = NA)
    lines(x.yrs,apply(est_mcmc,2,quantile,na.rm=T)[1,])
    lines(x.yrs,apply(est_mcmc,2,quantile,na.rm=T)[2,])
    lines(x.yrs,apply(est_mcmc,2,quantile,na.rm=T)[3,],lwd=3)
    lines(x.yrs,apply(est_mcmc,2,quantile,na.rm=T)[4,])
    lines(x.yrs,apply(est_mcmc,2,quantile,na.rm=T)[5,])
    #abline(v=c(1996, 2002, 2008),lty=2,lwd=2)
    #abline(v=c(2016),lty=2,lwd=2,col='blue')
    mtext('Year',outer=T,side=1,cex=1.4)
    mtext(paste(path.in,config.nm[i.cfg],sep=''),outer=T,side=3,cex=1.4)
   } #end plot
  } #i.est
  } #i.area
  } #i.cfg
graphics.off()
