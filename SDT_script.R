library("RTMB")

#suggestions for changes commented in below.

# could add to both functions below by incorporating option to sample from distributions
#as suggested by HS
# note would have to constrain criteria to be strictly increasing if these
#are not just fixed values.
#also can add sigma and option for evsd (sig=1) or uvsd for each of the
#probabiltiy and generating functions. would need to correct others to also
#have sigma with its own sd. 
#also might be nice to have this more flexible for any number of criteria
#if that works can then experiment with different parameterizations
gen_pars<-function(){
  list(dprime=.9,sd_ldprime=.2,crit=c(crit1=-1.7,crit2=0,crit3=.9))
}
init_pars<-function(){
  list(dprime=.1,sd_ldprime=.4,crit=c(crit1=-1.5,crit2=.2,crit3=.75))
}

optim_pars<-function(par,nT){
  #here constrained on log to be positive in optimizer 
  #and reverted to right scale below, but if randomly drawn values
  #then need to check for zeros
  mu_ldprime<-log(par$dprime)
  list(mu_ldprime=mu_ldprime,lsd_ldprime=log(par$sd_ldprime),
       crit1=unname(par$crit[1]),lgap=log(diff(par$crit)),
       ldprime_i=rep(mu_ldprime,nT))
}

evsd_prob<-function(dprime,crit){
  p<-c(pnorm(crit[1],dprime,1),diff(pnorm(crit,dprime,1)),
       1-pnorm(crit[length(crit)],dprime,1))
  p/sum(p)
}

#could change this to  loop through people n and trial n
sim_data<-function(sim,nT=125,nC=150,gen_p=gen_pars()){
  seed<-sample.int(.Machine$integer.max,1);set.seed(seed)
  n_sig<-n_noz<-nC/2
  n_bin<-length(gen_p$crit)+1
  #would need additional draws for sigma
  ldprime_i<-rnorm(nT,log(gen_p$dprime),gen_p$sd_ldprime)
  #note that this noise is the same for all people because
  #they have the same criteria. would be good to make this more 
  #flexible
  noz_prob<-evsd_prob(0,gen_p$crit)
  #convert back here to right scale
  #again would need to modify below for uvsd
  sig_cts<-t(vapply(exp(ldprime_i),function(x)
    rmultinom(1,n_sig,evsd_prob(x,gen_p$crit))[,1],integer(n_bin)))
  noz_cts<-t(rmultinom(nT,n_noz,noz_prob))
  colnames(sig_cts)<-colnames(noz_cts)<-paste0("bin",seq_len(n_bin))
  #save for checks and comparisons
  list(sim=sim,seed=seed,nT=nT,nC=nC,n_sig=n_sig,n_noz=n_noz,
       sig_cts=sig_cts,noz_cts=noz_cts,ldprime_i=ldprime_i)
}

rtmb_obj<-function(sim_dat,start=init_pars()){
  dat_list<-list(sig_cts=sim_dat$sig_cts,noz_cts=sim_dat$noz_cts)
  par_list<-optim_pars(start,nrow(sim_dat$sig_cts))
  nll<-function(pars){
    getAll(dat_list,pars,warn=FALSE)
    r_dprime<-exp(mu_ldprime); r_sd_ldprime<-exp(lsd_ldprime); r_crit<-crit1+c(0,cumsum(exp(lgap)))
    noz_prob<-evsd_prob(0,r_crit)
    nllout<- -sum(dnorm(ldprime_i,mu_ldprime,r_sd_ldprime,log=TRUE))
    for(i in seq_len(nrow(sig_cts))){
      sig_prob<-evsd_prob(exp(ldprime_i[i]),r_crit)
      nllout<-nllout-sum(sig_cts[i,]*log(sig_prob+1e-12))
      nllout<-nllout-sum(noz_cts[i,]*log(noz_prob+1e-12))
    }
    ADREPORT(r_dprime);ADREPORT(r_sd_ldprime);ADREPORT(r_crit)
    nllout
  }
  MakeADFun(nll,par_list,random="ldprime_i")
}

group_rec<-function(sdr,gen_p=gen_pars()){
  savtab<-as.data.frame(summary(sdr,"report"))
  gen_vals<-c(dprime=gen_p$dprime,sd_ldprime=gen_p$sd_ldprime,gen_p$crit)
  out<-data.frame(par=names(gen_vals),gen=unname(gen_vals),estimate=savtab[[1]],se=savtab[[2]])
  out$bias<-out$estimate-out$gen
  out$ci_low<-out$estimate-1.96*out$se
  out$ci_high<-out$estimate+1.96*out$se
  out$covered<-out$ci_low<=out$gen&out$gen<=out$ci_high
  out
}

indv_rec<-function(obj,opt,sim_dat){
  obj$fn(opt$par);
  indl<-obj$env$parList(obj$env$last.par.best)
  data.frame(sim=sim_dat$sim,
             ind=seq_along(indl$ldprime_i),
             ldprime_gen=sim_dat$ldprime_i,ldprime_est=as.numeric(indl$ldprime_i),
             dprime_gen=exp(sim_dat$ldprime_i),dprime_est=exp(as.numeric(indl$ldprime_i)))
}

dat<-sim_data(sim=1);
obj<-rtmb_obj(dat)
start_nll<-obj$fn(obj$par);
start_grad<-obj$gr(obj$par)
opt<-nlminb(obj$par,obj$fn,obj$gr);
sdr<-sdreport(obj,par.fixed=opt$par)

fit_check<-data.frame(
  code=opt$convergence,
  msg=opt$message,
  start_nll=start_nll,
  final_nll=opt$objective,
  max_start_grad=max(abs(start_grad)),
  pdHess=isTRUE(sdr$pdHess))

recovered<-group_rec(sdr,gen_pars())
ind_rec<-indv_rec(obj,opt,dat)

fit_check;recovered
head(ind_rec)
cor(ind_rec$dprime_gen,ind_rec$dprime_est)
plot(ind_rec$dprime_gen,ind_rec$dprime_est,xlab="generative ind dprime",
     ylab="recovered ind dprime");
abline(0,1)
