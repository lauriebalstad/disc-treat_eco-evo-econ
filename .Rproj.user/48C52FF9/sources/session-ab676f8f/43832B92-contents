library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(viridis)

#----figure 2: economic outcomes----
# using base parameters

# set up & make function
{tr.vec <- c(0.35, 0.55, 0.75, 0.95) # treatment effect --> might be interesting to also do with wild fish to see how threshold & wild fish interact
  th.vec <- c(0.002, 1:18) # threshold
  # lossesF1_mat <- matrix(NA, nrow = length(th.vec), ncol = length(tr.vec))
  treat_num <- matrix(NA, nrow = length(th.vec), ncol = length(tr.vec))
  lice_burd <- matrix(NA, nrow = length(th.vec), ncol = length(tr.vec))
  
  State = NA # dummy parameter for running?
  
  econ <- function() {
    for (i in 1:length(tr.vec)) {
      
      t_s.0 = tr.vec[i]
      t_r.0 = t_s.0*0.25 # make 25% of t_s
      
      for (j in 1:length(th.vec)) {
        
        th.0 = th.vec[j]
        
        # redefine root & event with correct threshold and effecitivity
        rootfunc <- function(Time, y, Pars) {return (y[1] + y[2] - th.0*F_f.0)} # threshold hit
        eventfunc <- function(Time, y, Pars) {
          y[1] <- (1-t_r.0)*y[1]
          y[2] <- (1-t_s.0)*y[2]
          return(y)
        } # at the root, do this event
        
        Pars = c(l_s = l_s.0, c_s = c_s.0, l_r = l_r.0, c_r = c_r.0, 
                 u_s = u_s.0, u_r = u_r.0, m = m.0, mb = mb.0,
                 p_ff = p_ff.0, p_ll = p_ll.0, p_fl = p_fl.0, p_lf = p_lf.0, p_ww = p_ww.0, 
                 B_f = B_f.0, B_l = B_l.0, B_w = B_w.0,
                 Tr = Tr.0, t_r = t_r.0, t_s = t_s.0, b = b.0, th = th.0, h = h.0,
                 F_f = F_f.0, F_w = F_w.0, 
                 r = r.0, v = v.0, X = X.0, Y = Y.0, u_f = u_f.0, sig = sig.0,
                 M_out = M_out.0, M_in = M_in.0)
        
        State = c(L_rF = 0, L_sF = 0.001, L_rL = 0, L_sL = 0, F_l = r.0, L_rW = 0, L_sW = F_w.0)
        out <- ode(Time, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
        # then take output from out and use to start new simulation state
        State = out[dim(out)[1], 2:8]; State["L_sF"] = State["L_sF"]-0.001; State["L_rF"] = 0.001 # can just replace a suseptable louse with a resistant louse
        out <- ode(Time, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE, maxroot = 6000))
        
        out_df <- as.data.frame(out)
        treat_time <- attributes(out)$troot
        treat_apply <- length(attributes(out)$troot)
        
        treat_num[j,i] <- treat_apply
        lice_burd[j,i] <- sum(out_df$L_rF + out_df$L_sF) # bc this is the sum of two lists
        
      }
    }
    
    return(list(treat_num, lice_burd))
    
  }}

# run econ function
{
econ_vals <- econ()
}
# manipulate dataframes
{# melt & join outcomes
treatment_apps <- melt(econ_vals[[1]])
lice_burd <- melt(econ_vals[[2]])
for (i in 1:dim(treatment_apps)[1]) {
  treatment_apps$Treat_Int[i] <- tr.vec[treatment_apps$Var2[i]] # rename based on vectors
  treatment_apps$Threshold[i] <- th.vec[treatment_apps$Var1[i]] # rename based on vectors
  lice_burd$Treat_Int[i] <- tr.vec[lice_burd$Var2[i]] # rename based on vectors
  lice_burd$Threshold[i] <- th.vec[lice_burd$Var1[i]] # rename based on vectors
}
# cleaning
lice_burd$tot_burd <- lice_burd$value
treatment_apps$tot_treat <- treatment_apps$value
lice_burd <- lice_burd[, 4:6]
treatment_apps <- treatment_apps[, 4:6]
# check burden 
lice_burd$perfish <- lice_burd$tot_burd/F_f.0/length(Time) # recall summed across all time steps
# join data
econ_dat <- merge(treatment_apps, lice_burd, by = c("Treat_Int", "Threshold"))

# add col for low/high cost
# goal is that cost is ~0.2*4.5*F_f.0 = 5400000
econ_dat$'high treatment cost' <- 2e-4*F_f.0*econ_dat$tot_treat + 5e-6*econ_dat$tot_burd
econ_dat$'low treatment cost' <- 1e-4*F_f.0*econ_dat$tot_treat + 5e-6*econ_dat$tot_burd
# convert data
econ_long <- econ_dat %>% pivot_longer(cols = c('high treatment cost', 'low treatment cost'))
}

# make figures
{no_treat_loss <- econ_long$value[which(econ_long$Threshold == max(econ_long$Threshold))[1]]
fig2 <- ggplot(econ_long, aes(Threshold, value, col = as.factor(Treat_Int))) + 
  geom_hline(aes(yintercept = no_treat_loss), lty = 2, col = "gray") + 
  geom_point(size = 1.4) + 
  facet_wrap(~name) +
  scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
  labs(x = "Treatment threshold", y = expression("Economic losses ("~10^{6}~")"), col = "Treatment \nefficacy") +
  scale_y_continuous(breaks = c(1e6, 2e6, 3e6), labels = c(1, 2, 3))
}

fig2
png("plots/economic_outcomes.png",height=90,width=170,res=400,units='mm')
print(fig2)
dev.off()


#----figure 3: ecological outcomes----
