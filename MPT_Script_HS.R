
library("RTMB")

data(arnold2013, package = "TreeBUGS")
d.encoding <- subset(arnold2013, group == "encoding")

nsubj <- nrow(d.encoding)
parameters <- list(
  b_mu = 0,
  d2_mu = 0,
  D3_mu = 0,
  g_mu = 0 ,
  b_sd = 1,
  d2_sd = 1,
  D3_sd = 1,
  g_sd = 1 ,
  b_u = rep(0, nsubj),
  d2_u = rep(0, nsubj),
  D3_u = rep(0, nsubj),
  g_u = rep(0, nsubj)
)

dlist <- list(
  data1 = as.matrix(d.encoding[,c("EE", "EN", "EU")]),
  data2 = as.matrix(d.encoding[,c("NE", "NN", "NU")]),
  data3 = as.matrix(d.encoding[,c("UE", "UN", "UU")])
)
str(dlist)

fmulti <- function(parms) {
  getAll(dlist, parms, warn=FALSE)
  
  nll <- 0
  nll <- nll - sum(dnorm(b_u, mean=b_mu, sd=b_sd, log=TRUE))
  nll <- nll - sum(dnorm(d2_u, mean=d2_mu, sd=d2_sd, log=TRUE))
  nll <- nll - sum(dnorm(D3_u, mean=D3_mu, sd=D3_sd, log=TRUE))
  nll <- nll - sum(dnorm(g_u, mean=g_mu, sd=g_sd, log=TRUE))
  
  b <- pnorm(b_u)
  d2 <- pnorm(d2_u)
  D3 <- pnorm(D3_u)
  g <- pnorm(g_u)
  
  # prob1 <- matrix(0, nsubj, 3)
  # prob2 <- matrix(0, nsubj, 3)
  # prob3 <- matrix(0, nsubj, 3)
  
  prob1 <- D3*d2 + D3*(1-d2)*g + (1-D3)*b*g
  prob2 <- (1-D3)*(1-b)
  prob3 <- D3*(1-d2)*(1-g) + (1-D3)*b*(1-g)
  prob4 <- (1-D3)*b*g
  prob5 <- D3 + (1-D3)*(1-b)
  prob6 <- (1-D3)*b*(1-g)
  prob7 <- D3*(1-d2)*g + (1-D3)*b*g
  prob8 <- (1-D3)*(1-b)
  prob9 <- D3*d2 + D3*(1-d2)*(1-g) + (1-D3)*b*(1-g)
  
  for (i in 1:nrow(data1)) {
    # prob1[subj, 1] <- D1[subj]*d1[subj] + D1[subj]*(1-d1[subj])*a[subj] + (1-D1[subj])*b[subj]*g[subj]
    # prob1[subj, 2] <- (1-D1[subj])*(1-b[subj])
    # prob1[subj, 3] <- D1[subj]*(1-d1[subj])*(1-a[subj]) + (1-D1[subj])*b[subj]*(1-g[subj])
    # prob2[subj, 1] <- (1-D3[subj])*b[subj]*g[subj]
    # prob2[subj, 2] <- D3[subj] + (1-D3[subj])*(1-b[subj])
    # prob2[subj, 3] <- (1-D3[subj])*b[subj]*(1-g[subj])
    # prob3[subj, 1] <- D2[subj]*(1-d2[subj])*a[subj] + (1-D2[subj])*b[subj]*g[subj]
    # prob3[subj, 2] <- (1-D2[subj])*(1-b[subj])
    # prob3[subj, 3] <- D2[subj]*d2[subj] + D2[subj]*(1-d2[subj])*(1-a[subj]) + (1-D2[subj])*b[subj]*(1-g[subj])
    nll <- nll - dmultinom(x = data1[i,], prob = c(prob1[i],prob2[i],prob3[i]), log = TRUE)
    nll <- nll - dmultinom(x = data2[i,], prob = c(prob4[i],prob5[i],prob6[i]), log = TRUE)
    nll <- nll - dmultinom(x = data3[i,], prob = c(prob7[i],prob8[i],prob9[i]), log = TRUE)
  }
  nll
}

fmulti(parameters) ## just check if it works

obj <- MakeADFun(fmulti, parameters, random = paste0(c("b", "d2", "D3", "g"), "_u"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
#zapsmall(opt_multi$par)
sdr_multi <- sdreport(obj)
sdr_multi
#summary(sdr_multi)
round(sdr_multi$par.fixed, 2)
#  b_mu d2_mu D3_mu  g_mu  b_sd d2_sd D3_sd  g_sd
# -0.15  0.41 -0.68  0.26  0.48  0.00  0.10  0.59
### b, a/g, D3/D1 basically match below, d2/d1 is a bit off.


EQNfile <- system.file("MPTmodels/2htsm.eqn", package = "TreeBUGS")
fit <- TreeBUGS::traitMPT(EQNfile, d.encoding,
                          n.thin = 5,
                          restrictions = list("D1=D2=D3", "d1=d2", "a=g")
)
summary(fit)
# Mean/Median of latent-trait values (probit-scale) across individuals:
#                Mean    SD   2.5%    50%  97.5% Time-series SE n.eff  Rhat R_95%
# latent_mu_a   0.242 0.149 -0.059  0.242  0.531          0.007   523 1.002 1.004
# latent_mu_b  -0.150 0.119 -0.383 -0.149  0.086          0.005   555 1.034 1.115
# latent_mu_d1  0.608 0.517 -0.237  0.538  1.858          0.007  5596 1.001 1.003
# latent_mu_D1 -0.704 0.084 -0.879 -0.702 -0.546          0.001  3915 1.001 1.003
#
# Standard deviation of latent-trait values (probit scale) across individuals:
#                  Mean    SD  2.5%   50% 97.5% Time-series SE n.eff  Rhat R_95%
# latent_sigma_a  0.661 0.142 0.428 0.646 0.986          0.003  2220 1.005 1.016
# latent_sigma_b  0.531 0.102 0.364 0.519 0.760          0.002  2998 1.001 1.002
# latent_sigma_d1 1.010 1.169 0.027 0.626 4.358          0.027  1924 1.003 1.005
# latent_sigma_D1 0.152 0.100 0.008 0.138 0.382          0.003  1499 1.001 1.001