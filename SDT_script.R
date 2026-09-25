library("RTMB")

nT <- 125
nC <- 150
n_sig <- floor(nC/2)
n_noz <- nC - n_sig

gen_pars <- list(mu=.9, sigma = 1.2, sd_mu=.2, sd_sigma=.2, rho_log = .4,
                 crit=c(crit1=-1.7,crit2=0,crit3=.9), sd_crit = .25)

init_pars <- list(mu=.1, sigma = 1, sd_lmu=.4, sd_lsigma = .1, rho_log = 0,
                  crit=c(crit1=-1.5,crit2=.2,crit3=.75), sd_crit = .1)

#corr version
#gen_pars <- list(mu=.9, sigma = 1.2, sd_mu=.2, sd_sigma=.2, 
#                 crit=c(crit1=-1.7,crit2=0,crit3=.9), sd_crit = .25)

#init_pars <- list(mu=.1, sigma = 1, sd_lmu=.4, sd_lsigma = .1, 
#                  crit=c(crit1=-1.5,crit2=.2,crit3=.75), sd_crit = .1)

#cor_mat <- matrix(c(1, gen_pars$rho_log, gen_pars$rho_log, 1))
#sd_mat <- diag(c(gen_pars$sd_mu,gen_pars$sd_sigma))
#cov_mat <- sd_mat%*%cor_mat%*%sd_mat
#add_corr <- matrix(rnorm(nT*2),nT,2)%*%chol(cov_mat)

#lmu_pi <- log(gen_pars$mu) + add_corr[,1]
#lsigma_pi <- log(gen_pars$sigma) + add_corr[,1]

sdt_prob <- function(mu,sigma,crit){
    cum_p <- pnorm(crit,mu,sigma)
    p <- c(cum_p[1], diff(cum_p),pnorm((mu-crit[3])/sigma))
    p/sum(p)}

# participant effects
lmu_pi <- rnorm(nT, log(gen_pars$mu), gen$sd_lmu)
lsigma_pi <- rnorm(nT, log(gen_pars$sigma), gen$sd_lsigma)
c_diff_pi <- rnorm(nT, 0, gen$sd_crit)

mu_pi <- exp(lmu_pi)
sigma_pi <- exp(lsigma_pi)
crit_pi <- cbind( gen_par$crit[1] + c_diff_pi, gen_par$crit[2] + c_diff_pi, gen_par$crit[3] + c_diff_pi,)

sig_cts <- noz_cts <- matrix(0,nT,4)

#sample freqs
for (i in seq_len(nT)) 
  { temp_cnts_sig <- rmultinom(1, n_sig, sdt_prob(mu_pi[i],sigma_pi[i],crit_pi[i,]))
        sig_cts[i,] <- temp_cnts_sig [,1] 
    temp_cnts_noz <- rmultinom(1, n_sig, sdt_prob(0,1,crit_pi[i,]))
        noz_cts[i,] <- temp_cnts_noz [,1] 
   }

dat <- list(sig_cts, noz_cts)

pars <- list(m_lmu = log(init_pars$mu),m_lsigma=log(init_pars$sigma), log_sd_lmu=log(init_pars$sd_lmu), log_sd_lsigma=log(init_pars$sd_lsigma),
             crit1 = init_pars$crit[1], crit_diff = log(diff(init_pars%crit)), log_sd_crit = log(init_pars$sd_crit),
             lmu_pi = rep(log(init_pars$mu), nT), lsigma_pi = rep(log(init_pars$sigma),nT), crit_diff_pi = rep(0,nT)) 
             
             #,rho = 0)

nll <- function(pars) { 
    getAll(dat, pars, warn = FALSE)
    sd_mu <- exp(log_sd_lmu)
    sd_sigma <- exp(log_sd_lsigma)
    sd_crit <- exp(log_sd_crit)
    crit <- crit1 + c(0, cumsum(exp(crit_diff)))
  #  conv_rho <- rho     possible to restrict? use bounds
  #  Rmat <- matrix(c(1,cov_rho,cov_rho,1),2,2)
  #  Dmat <- diag(c(sd_mu,sd_sigma))
  #  CovMat <- Dmat%*%Rmat%*%Dmat
  # nllout <- -sum(RTMB::dmvnorm(cbind(lmu_pi, m_lmu), mu = c(m_lmu, m_lsigma), Sigma = CovMat, log = TRUE))
  
    nllout <- -sum(dnorm(lmu_pi, m_lmu, sd_mu, log = TRUE))
    nllout <- nllout -sum(dnorm(lsigma_pi, m_lsigma, sd_sigma, log = TRUE))
    nllout <- nllout -sum(dnorm(crit_diff_pi, 0, sd_crit, log = TRUE))
    for (i in seq_len(nrow(sig_cts))) {
        crit_ind <- crit + crit_diff_pi[i]
        mu_ind <- exp(lmu_pi[i])
        sigma_ind <- exp(lsigma_pi[i])
        sig_prob <- sdt_prob(mu_ind, sigma_ind, crit_ind)
        noz_prob <- sdt_prob(0, 1, crit_ind)
        sig_prob <- sig_prob + 1e-15
        noz_prob <- noz_prob + 1e-15
        sig_loglik <- sum(sig_cts[i,]*log(sig_prob))
        noz_loglik <- sum(noz_cts[i,]*log(noz_prob))
        nllout <- nllout - sig_loglik - noz_loglik
      }
    med_mu <- exp(m_lmu) 
    ADREPORT(med_mu)
    med_sigma <- exp(m_lsigma)
    ADREPORT(med_sigma)
    ADREPORT(sd_mu)
    ADREPORT(sd_sigma)
    ADREPORT(crit)
    ADREPORT(sd_crit)
  # ADREPORT(conv_rho)
 nllout
}

rtmb_obj<-MakeADFun(nll,pars,random=c("lmu_pi","lsigma_pi","crit_diff_pi"))

opt <- nlminb(rmtb_obj$par, rtmb_obj$fn, rtmb_obj$gr)

stdr <- sdreport(obj,par.fixed = opt$par)
print(summary(stdr,"report"))

print(opt$message)
print(opt$convergence)
print(stdr$pdHess)
