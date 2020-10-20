
#################################################################
# R Code:                                                       #
# Smoothing corrections for improving sample size recalculation # 
# rules in adaptive group sequential study designs              #
# (Carolin Herrmann, Geraldine Rauch)                           #
#################################################################

################################################################################
# To start the simulations, please:
# - specify your specific path to the program in Main.R
# - determine the design (first and initial incremental second stage sample size 
#   per group, global one-sided significance level (one sided), locally adjusted 
#   significance levels, futility stopping bound, factor for maximal sample size) 
#   in Main.R
# - decide on the number of simulations n_sim in Main.R
# - define the true underlying standardized effect size for the evaluation in 
#   Main.R
################################################################################


library(lattice)
library(RColorBrewer)
library(latticeExtra)
library(xtable)

#-----------------------------------------------#
# Define your path to the program:
#-----------------------------------------------#
setwd("C:\\Users\\cherrman\\Desktop\\MIM_Revision\\2020_10_19_Acceptance\\Code_Upload")

#-----------------------------------------------------------------------------#
# Determine the study design with
# first stage sample size per group n_1, 
# adjusted local significance level after first stage alpha_1,
# adjusted local significance level after second stage alpha_1_2,
# initially incremental second stage sample size per group n_2, 
# futility stop bound alpha_0,
# factor for maximal sample size f with n_max=f*n_1,
# global significance level (one sided) alpha_glob 
#-----------------------------------------------------------------------------#
design <- list(
  n_1 = 50,
  alpha_1 = 0.002638, # (O'Brien Fleming: 0.002638; Pocock: 0.01476; Wang Tsiatis with 0.25: 0.007766)
  alpha_1_2 = 0.024270, # (O'Brien Fleming: 0.024270; Pocock: 0.01476; Wang Tsiatis with 0.25: 0.020938)
  n_2 = 50,
  alpha_0 = 0.5,
  f = 4,
  alpha_glob = 0.025
)

#---------------------------------------------#
# Specify the number of simulations n_sim
#---------------------------------------------#
n_sim <- 10000

source("cond_power.r") 
source("cond_score.r") 

# Sample size recalculation according to "restricted observed conditional
# power approach" without smoothing and with five different smoothing options:
source("restrOCP_design.r") 
source("restrOCP_design_linear.r") 
source("restrOCP_design_threesteps.r")
source("restrOCP_design_sigmoid.r")
source("restrOCP_design_concave.r")
source("restrOCP_design_convex.r")


#-------------------------------------------------------------------------------
# Plot the sample size as a function of t1 and delta
#-------------------------------------------------------------------------------
set.seed(5)
t1 <- seq(-1,3,0.01)
delta <- t1*sqrt(2/design$n_1)
cincr <- (qnorm(1-design$alpha_1_2) * sqrt(design$n_1 + design$n_2)/sqrt(design$n_2) - qnorm(0.4)) /
  (sqrt(design$n_1)/sqrt(design$n_2) + sqrt((design$f*design$n_1 - design$n_1) / design$n_1))
cdecr <- (qnorm(1-design$alpha_1_2) * sqrt(design$n_1 + design$n_2)/sqrt(design$n_2) - qnorm(0.2)) /
  (sqrt(design$n_1)/sqrt(design$n_2) + sqrt((design$f*design$n_1 - design$n_1) / design$n_1))

par(mar = c(5.1, 4.1, 4.1, 2.1), cex.axis = 1, cex.lab = 1, mgp = c(3,1,0))
plot(restrOCP_design(design, delta/sqrt(2/design$n_1), delta)$n ~ delta, 
     xlim = c(-1*sqrt(2/design$n_1), 3*sqrt(2/design$n_1)), ylim = c(0,300), 
     type = 'l', col = "blue", lwd = 2, lty = 1, xlab = NA, ylab = NA, axes = F)
axis(side = 3)
mtext(side = 3, line = 3, 'Observed Interim Effect')
mtext(side = 1, line = 3, 'Observed Value of the Interim Test Statistic')
par(new = TRUE)
plot(restrOCP_design_linear(design, t1, delta)$n ~ t1, xlim = c(-1,3), ylim = c(0,300),
     type = 'l', col = "green2", lwd = 2, lty = 2, xlab = NA, ylab = "Total Sample Size per Group")
par(new = TRUE)
plot(restrOCP_design_threesteps(design, t1, delta)$n ~ t1, xlim = c(-1,3), ylim = c(0,300),
     type = 'l', col = "purple", lwd = 2, lty = 4, xlab = NA, ylab = NA)
par(new = TRUE)
plot(restrOCP_design_sigmoid(design, t1, delta, 10)$n ~ t1, xlim = c(-1,3), ylim = c(0,300),
     type = 'l', col = "magenta", lwd = 2, lty = 4, xlab = NA, ylab = NA)
par(new = TRUE)
plot(restrOCP_design_concave(design, t1, delta)$n ~ t1, xlim = c(-1,3), ylim = c(0,300),
     type = 'l', col = "orange2", lwd = 2, lty = 3, xlab = NA, ylab = NA)
par(new = TRUE)
plot(restrOCP_design_convex(design, t1, delta)$n ~ t1, xlim = c(-1,3), ylim = c(0,300),
     type = 'l', col = "grey1", lwd = 2, lty = 3, xlab = NA, ylab = NA)
par(new = TRUE)
abline(v = qnorm(1-design$alpha_0), col="black", lty = 1)
text(qnorm(1-design$alpha_0)-0.06, 5, expression(c[fut]))
text(cincr-0.1, 5, expression(c[incr]))
text(cdecr+0.11, 5, expression(c[decr]))
text(qnorm(1-design$alpha_1)+0.06, 5, expression(c[eff]))
abline(v = qnorm(1-design$alpha_1), col="black", lty = 1)
abline(v = cincr, col="black", lty = 3)
abline(v = cdecr, col="black", lty = 3)
par(new = TRUE)
legend(-1.15, 300, c("Without", "Linear", "Stepwise", "Sigmoid", "Concave", "Convex"),
       lty = c(1, 2, 4, 4, 3, 3), lwd = c(2, 2, 2, 2, 2, 2), bty = "n", 
       col = c("blue", "green2", "purple", "magenta", "orange2", "grey1"))


###################################################################
# Initiate Performance Simulation
###################################################################

#------------------------------------------------------------
# Specify the underlying treatment effect Delta (e.g. 0.4)
#------------------------------------------------------------
delta <- 0.4

#------------------------------------------------------------
# Create values of test statistic
#------------------------------------------------------------
set.seed(140)
s1 <- rnorm(n_sim, delta*sqrt(design$n_1/2), 1)

#--------------------------------------------------------------------
# Determine required fixed sample size
# and define the target values for sample size and conditional power
#--------------------------------------------------------------------
if(delta != 0){
  n_req <- 1 * ((qnorm(0.975) + qnorm(0.8)) / delta)^2
  pow <- power.t.test(n = n_req, delta = delta, sd = 1, sig.level = 0.025,
                    power = NULL,
                    type = c("two.sample"),
                    alternative = "one.sided")$power
  while(pow < 0.8){
    n_req <- n_req + 1
    pow <- power.t.test(n = n_req, delta = delta, sd = 1, sig.level = 0.025,
                       power = NULL,
                       type = c("two.sample"),
                       alternative = "one.sided")$power
  }
  n_req <- n_req
  ref_cp <- 0.8
}
if(delta == 0){
  n_req <- design$n_1
  ref_cp <- 0.025
}

if(n_req > design$f * design$n_1){
  n_req <- design$n_1
  ref_cp <- 0.025
}

#---------------------------------------------------------------------------
# Restrict the simulated data to the recalculation area [c_fut; c_eff]
#---------------------------------------------------------------------------
t1 <- s1[s1 < qnorm(1 - design$alpha_1) & s1 >= qnorm(1 - design$alpha_0)]

#------------------------------------------------------------------------
# Store the recalculated total sample sizes "n" per design in one matrix
#------------------------------------------------------------------------
n <- matrix(nrow = 6, ncol = length(t1))
rownames(n) <- c("restrOCP", "restrOCPlinear", "restrOCPsteps", "restrOCPsigmoid", "restrOCPconcave", "restrOCPconvex")
n[1,] <- restrOCP_design(design, t1, delta)$n
n[2,] <- restrOCP_design_linear(design, t1, delta)$n
n[3,] <- restrOCP_design_threesteps(design, t1, delta)$n
n[4,] <- restrOCP_design_sigmoid(design, t1, delta, steep = 10)$n
n[5,] <- restrOCP_design_concave(design, t1, delta)$n
n[6,] <- restrOCP_design_convex(design, t1, delta)$n
#--------------------------------------------------------------------------------
# Store the conditional power "cp" for the observed delta at the interim analysis
#--------------------------------------------------------------------------------
par(mfrow = c(1,1))
cp <- matrix(nrow = 6, ncol = length(t1))
rownames(cp) <- c("restrOCP", "restrOCPlinear", "restrOCPsteps", "restrOCPsigmoid", "restrOCPconcave", "restrOCPconvex")
cp[1,] <- restrOCP_design(design, t1, delta)$cp
cp[2,] <- restrOCP_design_linear(design, t1, delta)$cp
cp[3,] <- restrOCP_design_threesteps(design, t1, delta)$cp
cp[4,] <- restrOCP_design_sigmoid(design, t1, delta, steep = 10)$cp
cp[5,] <- restrOCP_design_concave(design, t1, delta)$cp
cp[6,] <- restrOCP_design_convex(design, t1, delta)$cp
#-----------------------------------------------------------
# Store the conditional power "cp_true" for the true delta 
#-----------------------------------------------------------
par(mfrow = c(1,1))
cp_true <- matrix(nrow = 6, ncol = length(t1))
rownames(cp_true) <- c("restrOCP", "restrOCPlinear", "restrOCPsteps", "restrOCPsigmoid", "restrOCPconcave", "restrOCPconvex")
cp_true[1,] <- restrOCP_design(design, t1, delta)$cp_true
cp_true[2,] <- restrOCP_design_linear(design, t1, delta)$cp_true
cp_true[3,] <- restrOCP_design_threesteps(design, t1, delta)$cp_true
cp_true[4,] <- restrOCP_design_sigmoid(design, t1, delta, steep = 10)$cp_true
cp_true[5,] <- restrOCP_design_concave(design, t1, delta)$cp_true
cp_true[6,] <- restrOCP_design_convex(design, t1, delta)$cp_true

#----------------------------------------------------------------------------
# output for conditional performance out_cond:
# expected sample size per group when entering the recalculation area 
# (RA) "Mean n", variance of sample size per group when entering the 
# RA "Var n", expected conditional power when entering the RA "Mean cp", 
# variance of conditional power when entering the RA "Var cp", point-
# wise conditional score by Herrmann et al. (2020) "cond score"
#----------------------------------------------------------------------------
nscore <- cond_score(design, cp, ref_cp, n, n_req, n_sim)
out_cond <- matrix(nrow = 6, ncol = 5)
rownames(out_cond) <- c("Without", "Linear", "Stepwise", "Sigmoid", "Concave", "Convex")
colnames(out_cond) <- c("Mean n", "Var n", "Mean cp", "Var cp", 
                        "cond score")  
for(i in 1:6){
  out_cond[i,] <- round(c(nscore$mean_n[i], nscore$var_n[i], 
                          nscore$mean_cp[i], nscore$var_cp[i], 
                          nscore$score_cond[i]), digits = 3)
}
print(out_cond)
