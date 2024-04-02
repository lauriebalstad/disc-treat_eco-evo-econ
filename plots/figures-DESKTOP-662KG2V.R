library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(viridis)
library(doParallel)
library(cowplot)

# economic outcomes: watch number of roots.... maxroot = 1000? then, filter to only last month?? 
# play with parameters: thinking going "standardized" to give sense of patterns, shapes, trends more than accuracy?
# need to standardize outputs too
# and then watch minimum trigger, too low right now

#----figure 2: economic outcomes----
# using base parameters

# set up & make function
{
tr.vec <- seq(from = 0.2, to = 0.9, length.out = 5) # treatment effect --> might be interesting to also do with wild fish to see how threshold & wild fish interact
th.vec <- seq(0.01, 3.01, length.out = 16) # threshold
# tr.vec <- c(0.35, 0.55, 0.75, 0.95) # treatment effect --> might be interesting to also do with wild fish to see how threshold & wild fish interact
# th.vec <- c(0.0015, 1:18) # threshold
# lossesF1_mat <- matrix(NA, nrow = length(th.vec), ncol = length(tr.vec))
treat_num <- matrix(NA, nrow = length(th.vec), ncol = length(tr.vec))
lice_burd <- matrix(NA, nrow = length(th.vec), ncol = length(tr.vec))

State = NA # dummy parameter for running?

econ <- function() {
  for (i in 1:length(tr.vec)) {
    
    t_s.0 = tr.vec[i]
    t_r.0 = t_s.0*0.5 # make 25% of t_s
    
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
               t_r = t_r.0, t_s = t_s.0, th = th.0, h = h.0,
               F_f = F_f.0, F_w = F_w.0, 
               r = r.0, v = v.0, X = X.0, Y = Y.0, u_f = u_f.0, sig = sig.0,
               M_out = M_out.0, M_in = M_in.0)
      
      State = c(L_rF = 0, L_sF = 0.001, L_rL = 0, L_sL = 0, F_l = r.0, L_rW = 0, L_sW = F_w.0)
      out <- ode(Time, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
      # then take output from out and use to start new simulation state
      State = out[dim(out)[1], 2:8]; State["L_sF"] = State["L_sF"]-0.001; State["L_rF"] = 0.001 # can just replace a suseptable louse with a resistant louse
      out <- ode(Time, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE, maxroot = 10000))
      
      # just last month of things
      out_df <- as.data.frame(out) %>% filter(time > Time-31)
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
# manipulate dataframes to add economics columns
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
  econ_dat$'high treatment cost' <- 2.5e-4*F_f.0*econ_dat$tot_treat + 5e-6*econ_dat$tot_burd
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

# using base parameters

# set up & make function
{
tr.vec <- seq(from = 0.2, to = 0.9, length.out = 5) # treatment effect --> might be interesting to also do with wild fish to see how threshold & wild fish interact
th.vec <- seq(0.01, 3.01, length.out = 16) # threshold
# tr.vec <- c(0.35, 0.55, 0.75, 0.95) # treatment effect --> might be interesting to also do with wild fish to see how threshold & wild fish interact
# th.vec <- c(0.0015, 1:18) # threshold
# lossesF1_mat <- matrix(NA, nrow = length(th.vec), ncol = length(tr.vec))
link_fish <- matrix(NA, nrow = length(th.vec), ncol = length(tr.vec))

State = NA # dummy parameter for running?

linkpop <- function() {
  for (i in 1:length(tr.vec)) {
    
    t_s.0 = tr.vec[i]
    t_r.0 = t_s.0*0.5 # make 25% of t_s
    
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
              t_r = t_r.0, t_s = t_s.0, th = th.0, h = h.0,
               F_f = F_f.0, F_w = F_w.0, 
               r = r.0, v = v.0, X = X.0, Y = Y.0, u_f = u_f.0, sig = sig.0,
               M_out = M_out.0, M_in = M_in.0)
      
      State = c(L_rF = 0, L_sF = 0.001, L_rL = 0, L_sL = 0, F_l = r.0, L_rW = 0, L_sW = F_w.0)
      out <- ode(Time, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
      # then take output from out and use to start new simulation state
      State = out[dim(out)[1], 2:8]; State["L_sF"] = State["L_sF"]-0.001; State["L_rF"] = 0.001 # can just replace a suseptable louse with a resistant louse
      out <- ode(Time, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE, maxroot = 6000))
      
      out_df <- as.data.frame(out)
      # get the last 30 days
      out_last_month <- out_df[(length(Time)-30):length(Time), ]
      
      link_fish[j,i] <- mean(out_last_month$F_l)
        
    }
  }
  
  return(link_fish)
  
}}

# run econ function
{
  link_vals <- linkpop()
}
# manipulate dataframes
{# melt & join outcomes
  link_size <- melt(link_vals)
  for (i in 1:dim(link_size)[1]) {
    link_size$Treat_Int[i] <- tr.vec[link_size$Var2[i]] # rename based on vectors
    link_size$Threshold[i] <- th.vec[link_size$Var1[i]] # rename based on vectors
  }
}

# make figures
{fig3 <- ggplot(link_size, aes(Threshold, value, col = as.factor(Treat_Int))) + 
    geom_point(size = 1.4) + 
    scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
    labs(x = "Treatment threshold", y = expression("Link Population ("~10^{3}~")"), col = "Treatment \nefficacy") +
    scale_y_continuous(breaks = c(2e3, 3e3, 4e3, 5e3), labels = c(2, 3, 4, 5))
}

fig3
png("plots/ecological_outcomes.png",height=90,width=85,res=400,units='mm')
print(fig3)
dev.off()

#----figure 4: evolution outcomes----
# could add a doparallel call... takes a while
# added a few additional efficacies to show range better?
# might need to modify base parameters so that HDRE is a real thing???

# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

# set up & make function
tr.vec <- seq(from = 0.2, to = 0.9, length.out = 5) # treatment effect --> might be interesting to also do with wild fish to see how threshold & wild fish interact
th.vec <- seq(0.01, 3.01, length.out = 16) # threshold
tr_th <- expand.grid(treat_eff = tr.vec, thresh_val = th.vec)
tr_th$connectivity <- rep(p_ff.0); tr_th$wildfish <- rep(F_w.0); tr_th$productivity <- rep(r.0); tr_th$type <- rep("base parameters")
tr_th_conn <- expand.grid(treat_eff = tr.vec, thresh_val = th.vec)
tr_th_conn$connectivity <- rep(p_ff.0-0.15); tr_th_conn$wildfish <- rep(F_w.0); tr_th_conn$productivity <- rep(r.0); tr_th_conn$type <- rep("increased connectivity")
tr_th_refg <- expand.grid(treat_eff = tr.vec, thresh_val = th.vec)
tr_th_refg$connectivity <- rep(p_ff.0); tr_th_refg$wildfish <- rep(F_w.0*3); tr_th_refg$productivity <- rep(r.0); tr_th_refg$type <- rep("increased refuge")
tr_th_prod <- expand.grid(treat_eff = tr.vec, thresh_val = th.vec)
tr_th_prod$connectivity <- rep(p_ff.0); tr_th_prod$wildfish <- rep(F_w.0); tr_th_prod$productivity <- rep(r.0*3); tr_th_prod$type <- rep("increased link productivity")
TR_TH <- rbind(tr_th, tr_th_conn, tr_th_refg, tr_th_prod)
# TR_TH <- rbind(tr_th)

# think about how to prep the foreach loop --> check how you did it in the evo rescue files?
base_pars = c(l_s = l_s.0, c_s = c_s.0, l_r = l_r.0, c_r = c_r.0, 
         u_s = u_s.0, u_r = u_r.0, m = m.0, mb = mb.0,
         p_ff = p_ff.0, p_ll = p_ll.0, p_fl = p_fl.0, p_lf = p_lf.0, p_ww = p_ww.0, 
         B_f = B_f.0, B_l = B_l.0, B_w = B_w.0,
         t_r = t_r.0, t_s = t_s.0, th = th.0, h = h.0,
         F_f = F_f.0, F_w = F_w.0, 
         r = r.0, v = v.0, X = X.0, Y = Y.0, u_f = u_f.0, sig = sig.0,
         M_out = M_out.0, M_in = M_in.0)

resist_percent_list <- foreach(i=1:dim(TR_TH)[1], .packages="deSolve") %dopar% {
    
      t_s.0 = TR_TH$treat_eff[i]
      t_r.0 = t_s.0*0.5 # make 25% of t_s
      p_ff.0 = TR_TH$connectivity[i] # low connectivity case (p_ff.0 close to 1) 
      p_ll.0 = p_ff.0 # lice in farm stay on farm...
      p_fl.0 = 1-p_ff.0 # ...or jump to wild (spillover)
      p_lf.0 = p_fl.0 # spillback rate
      p_ww.0 = 1 
      
      base_pars['t_s'] <- t_s.0; base_pars['t_r'] <- t_r.0
      base_pars['p_ff'] <- p_ff.0; base_pars['p_ll'] <- p_ll.0; base_pars['p_fl'] <- p_fl.0; base_pars['p_lf'] <- p_lf.0; base_pars['p_ww'] <- p_ww.0
      base_pars['r'] <- TR_TH$productivity[i]
      F_w.0 <- TR_TH$wildfish[i]
      
      th.0 = TR_TH$thresh_val[i]
      
      # redefine root & event with correct threshold and effecitivity
      rootfunc <- function(Time, y, base_pars) {return (y[1] + y[2] - th.0*F_f.0)} # threshold hit
      eventfunc <- function(Time, y, base_pars) {
        y[1] <- (1-t_r.0)*y[1]
        y[2] <- (1-t_s.0)*y[2]
        return(y)
      } # at the root, do this event
      
      # Pars['th.0']
      
      State = c(L_rF = 0, L_sF = 0.001, L_rL = 0, L_sL = 0, F_l = TR_TH$productivity[i], L_rW = 0, L_sW = F_w.0)
      out <- ode(Time, y = State, func = systemTREAT, parms = base_pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
      # then take output from out and use to start new simulation state
      State = out[dim(out)[1], 2:8]; State["L_sF"] = State["L_sF"]-0.001; State["L_rF"] = 0.001 # can just replace a suseptable louse with a resistant louse
      out <- ode(Time, y = State, func = systemTREAT, parms = base_pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE, maxroot = 6000))
      
      out_df <- as.data.frame(out)
      # get the last 30 days
      out_last_month <- out_df[(length(Time)-30):length(Time), ]
      
      resist_percent <- mean(out_last_month$L_rF/(out_last_month$L_rF + out_last_month$L_sF))
   
      }

# end parallel cluster
stopCluster(cl)

# manipulate dataframes
{
evol_dat <- TR_TH
evol_dat$value <- unlist(resist_percent_list)
}

# make figures
{fig4 <- ggplot(evol_dat, aes(thresh_val, value, col = as.factor(treat_eff))) + 
    geom_point(size = 1.4) + 
    facet_wrap(~type, nrow=2) +
    scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
    labs(x = "Treatment threshold", y = "Proportion resistant", col = "Treatment \nefficacy") 
}

fig4 # has some more efficacy options
png("plots/evolution_outcomes.png",height=170,width=170,res=400,units='mm')
print(fig4)
dev.off()

#----figure 5: frontiers----

# need high treatment cost & link population size
tmp <- merge(econ_dat[,c(1,2,6)], link_size[,c(3:5)], by = c("Treat_Int", "Threshold"))
colnames(tmp) <- c("treat_eff", "thresh_val", "econ", "ecol")
# need base parameter case
tmp2 <- evol_dat %>% filter(type == "increased connectivity") %>% select(c("treat_eff", "thresh_val", "value"))
colnames(tmp2) <- c("treat_eff", "thresh_val", "evol")
fronts <- merge(tmp, tmp2, by = c("treat_eff", "thresh_val"))

# need to add lines
ecol_econ <- ggplot(fronts, aes(ecol, -econ, col = as.factor(treat_eff))) + 
  geom_point(size = 1.4) + 
  scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
  labs(x = "Link population size", y = "Neg. economic costs", col = "Treatment \nefficacy") 
evol_econ <- ggplot(fronts, aes(-evol, -econ, col = as.factor(treat_eff))) + 
  geom_point(size = 1.4) + 
  scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
  labs(x = "Neg. prop. resist.", y = "Neg. economic costs", col = "Treatment \nefficacy") 
ecol_evol <- ggplot(fronts, aes(ecol, -evol, col = as.factor(treat_eff))) + 
  geom_point(size = 1.4) + 
  scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
  labs(x = "Link population size", y = "Neg. prop. resist.", col = "Treatment \nefficacy") 

ecol_econ
evol_econ
ecol_evol
