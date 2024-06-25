library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(viridis)
library(doParallel)
library(cowplot)

#----figure 2: economic outcomes----
# using base parameters

# set up & make function
{
  ep.vec <- seq(from = 0.3, to = 0.9, length.out = 3) # treatment effect --> might be interesting to also do with wild fish to see how threshold & wild fish interact
  M.vec <- seq(from = 0.01, to = 1.61, length.out = 9) # 13 -- threshold
  # ep.vec <- c(0.35, 0.55, 0.75, 0.95) # treatment effect --> might be interesting to also do with wild fish to see how threshold & wild fish interact
# M.vec <- c(0.0015, 1:18) # threshold
# lossesF1_mat <- matrix(NA, nrow = length(M.vec), ncol = length(ep.vec))
treat_num <- matrix(NA, nrow = length(M.vec), ncol = length(ep.vec))
lice_burd <- matrix(NA, nrow = length(M.vec), ncol = length(ep.vec))
link_pop <- matrix(NA, nrow = length(M.vec), ncol = length(ep.vec))

State = NA # dummy parameter for running?

outcomes_ecol_econ <- function(noLink) {
  
  if (noLink) {n_N1.0 <- 0} # if no link fish, remove reproduction
  
  for (i in 1:length(ep.vec)) {
    
    ep_s.0 = ep.vec[i]
    ep_r.0 = ep_s.0*0.05 # make 25% of ep_s
    
    for (j in 1:length(M.vec)) {
      
      M.0 = M.vec[j]
      # print(c(i, j))

      # redefine root & event with correct threshold and effecitivity
      rootfunc <- function(Time, y, Pars) {return (y[1] + y[2] - M.0*F_d.0)} # threshold hit
      eventfunc <- function(Time, y, Pars) {
        y[1] <- (1-ep_r.0)*y[1]
        y[2] <- (1-ep_s.0)*y[2]
        return(y)
      } # at the root, do this event
      
      Pars = base_pars
      Pars['ep_s'] <- ep_s.0
      Pars['ep_r'] <- ep_r.0
      
      State = c(L_rd = 0, L_sd = 0.001, L_rj = 0, L_sj = 0, F_j = n_N1.0, F_a = n_N1.0, L_rn = n_N1.0, L_sn = n_N1.0, F_n = n_N1.0, L_rw = 0, L_sw = F_w.0)
      out <- ode(Time, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
      # then take output from out and use to start new simulation state
      State = out[dim(out)[1], 2:2:dim(out)[2]]; State["L_sd"] = State["L_sd"]-0.001; State["L_rd"] = 0.001 # can just replace a suseptable louse with a resistant louse
      Time_invade = seq(0, 5*n-31, length.out = 5*n-30)
      out <- ode(Time_invade, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
      Time_outs = seq(0, 30, length.out = 31) # simulate just last month
      State = out[dim(out)[1], 2:dim(out)[2]]
      out <- ode(Time_outs, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE, maxroots = 20000))
      
      # just last month of things
      out_df <- as.data.frame(out) # %>% filter(time > (max(Time)-31))
      treat_time <- attributes(out)$troot
      treat_apply <- length(treat_time) # (which(treat_time > (max(Time)-31)))
      
      treat_num[j,i] <- treat_apply
      lice_burd[j,i] <- sum(out_df$L_rd + out_df$L_sd) # bc this is the sum of two lists
      link_pop[j,i] <- mean(out_df$F_j)
      
    }
  }
  
  return(list(treat_num, lice_burd, link_pop))
  
}}

# run econ function
{
  outcome_vals <- outcomes_ecol_econ(FALSE)
}

# manipulate dataframes to add economics columns -- for high/low...
{# melt & join outcomes
  treatment_apps <- melt(outcome_vals[[1]])
  lice_burd <- melt(outcome_vals[[2]])
  for (i in 1:dim(treatment_apps)[1]) {
    treatment_apps$Treat_Int[i] <- ep.vec[treatment_apps$Var2[i]] # rename based on vectors
    treatment_apps$Threshold[i] <- M.vec[treatment_apps$Var1[i]] # rename based on vectors
    lice_burd$Treat_Int[i] <- ep.vec[lice_burd$Var2[i]] # rename based on vectors
    lice_burd$Threshold[i] <- M.vec[lice_burd$Var1[i]] # rename based on vectors
  }
  # cleaning
  lice_burd$tot_burd <- lice_burd$value
  treatment_apps$tot_treat <- treatment_apps$value
  lice_burd <- lice_burd[, 4:6]
  treatment_apps <- treatment_apps[, 4:6]
  # check burden 
  lice_burd$perfish <- lice_burd$tot_burd/F_d.0/30 # length(Time) # recall summed across all time steps
  # join data
  econ_dat <- merge(treatment_apps, lice_burd, by = c("Treat_Int", "Threshold"))
  
  # add col for low/high cost
  # goal is that cost is ~0.3*4.5*F_d.0 = 8100000
  econ_dat$'High treatment cost' <- 0.075*F_d.0*econ_dat$tot_treat + .15*econ_dat$tot_burd
  econ_dat$'Low treatment cost' <- 0.025*F_d.0*econ_dat$tot_treat + .15*econ_dat$tot_burd # tot_burd
  max_econ <- max(c(econ_dat$'High treatment cost',econ_dat$'Low treatment cost'))
  econ_dat$'High treatment cost' <- econ_dat$'High treatment cost'/max_econ
  econ_dat$'Low treatment cost' <- econ_dat$'Low treatment cost'/max_econ
  # convert data
  econ_long <- econ_dat %>% pivot_longer(cols = c('High treatment cost', 'Low treatment cost'))
}

# make list of economic things
econ_outs <- list(base_pars, outcome_vals, econ_dat, econ_long)
# save(econ_outs, file = "simulated data/econ_outs.Rdata")

# make figures
{no_treat_loss <- econ_long$value[which(econ_long$Threshold == max(econ_long$Threshold))[1]]
  fig2 <- ggplot(econ_long, aes(Threshold, (value), col = as.factor(Treat_Int))) + 
    geom_hline(aes(yintercept = (no_treat_loss)), lty = 2, col = "gray") + 
    geom_point(size = 1.4) + # geom_line() +
    facet_wrap(~name, nrow = 2)  + #, scale = "free_y") +
    scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
    coord_cartesian(ylim = c(0, 1)) + 
    theme(legend.position = "bottom",
          legend.title = element_text(size = 10), 
          legend.text = element_text(size = 10), 
          axis.title = element_text(size = 10), 
          axis.text = element_text(size = 10), 
          strip.text = element_text(size = 10)) +
  labs(x = "Treatment threshold", y = expression("Economic losses (Standardized)"), col = "Treatment \nefficacy")
}

fig2

#----figure 3: ecological outcomes----

# using base parameters -- grab final data frame from before
# manipulate dataframes
{# melt & join outcomes
  link_size <- melt(outcome_vals[[3]])
  for (i in 1:dim(link_size)[1]) {
    link_size$Treat_Int[i] <- ep.vec[link_size$Var2[i]] # rename based on vectors
    link_size$Threshold[i] <- M.vec[link_size$Var1[i]] # rename based on vectors
  }
  link_size$standard <- link_size$value/max(link_size$value)
}

ecol_outs <- list(base_pars, outcome_vals, link_size)
# save(ecol_outs, file = "simulated data/ecol_outs.Rdata")

# make figures
{fig3 <- ggplot(link_size, aes(Threshold, standard, col = as.factor(Treat_Int))) + 
    geom_point(size = 1.4) + 
    scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
    theme(legend.position = c(0.85, 0.7),
          legend.title = element_text(size = 10), 
          legend.text = element_text(size = 10), 
          axis.title = element_text(size = 10), 
          axis.text = element_text(size = 10)) +
    labs(x = "Treatment threshold", y = "Juv. salmon population (Standardized)", col = "Treatment \nefficacy") 
}

fig3

#----figure 4: evolution outcomes----

# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

# set up & make function
{ep.vec <- seq(from = 0.3, to = 0.9, length.out = 3) # treatment effect --> might be interesting to also do with wild fish to see how threshold & wild fish interact
  M.vec <- seq(from = 0.01, to = 1.61, length.out = 9) # 13 -- threshold
  tr_th <- expand.grid(treat_eff = ep.vec, thresh_val = M.vec)
tr_th$connectivity <- rep(p_dd.0); tr_th$wildfish <- rep(F_w.0); tr_th$productivity <- rep(n_N1.0); tr_th$type <- rep("Base parameters")
tr_th_conn <- expand.grid(treat_eff = ep.vec, thresh_val = M.vec)
tr_th_conn$connectivity <- rep(p_dd.0-0.05); tr_th_conn$wildfish <- rep(F_w.0); tr_th_conn$productivity <- rep(n_N1.0); tr_th_conn$type <- rep("Increased connectivity")
tr_th_refg <- expand.grid(treat_eff = ep.vec, thresh_val = M.vec)
tr_th_refg$connectivity <- rep(p_dd.0); tr_th_refg$wildfish <- rep(F_w.0*2); tr_th_refg$productivity <- rep(n_N1.0); tr_th_refg$type <- rep("Increased refuge")
tr_th_prod <- expand.grid(treat_eff = ep.vec, thresh_val = M.vec)
tr_th_prod$connectivity <- rep(p_dd.0); tr_th_prod$wildfish <- rep(F_w.0); tr_th_prod$productivity <- rep(n_N1.0*1.25); tr_th_prod$type <- rep("Increased link productivity")
TR_TH <- rbind(tr_th, tr_th_conn, tr_th_refg, tr_th_prod)
# TR_TH <- rbind(tr_th, tr_th_refg)

# think about how to prep the foreach loop --> check how you did it in the evo rescue files?

resist_percent_list <- foreach(i=1:dim(TR_TH)[1], .packages="deSolve") %dopar% {
    
      ep_s.0 = TR_TH$treat_eff[i]
      ep_r.0 = ep_s.0*0.05 # make 75% of ep_s
      p_dd.0 = TR_TH$connectivity[i] # low connectivity case (p_dd.0 close to 1) 
      p_ll.0 = p_dd.0 # lice in farm stay on farm...
      p_dl.0 = 1-p_dd.0 # ...or jump to wild (spillover)
      p_ld.0 = p_dl.0 # spillback rate
      p_ww.0 = 1 
      
      Pars = base_pars
      Pars['ep_s'] <- ep_s.0; Pars['ep_r'] <- ep_r.0
      Pars['p_ff'] <- p_dd.0; Pars['p_ll'] <- p_ll.0; Pars['p_fl'] <- p_dl.0; Pars['p_lf'] <- p_ld.0; Pars['p_ww'] <- p_ww.0
      Pars['n_N1'] <- TR_TH$productivity[i]
      F_w.0 <- TR_TH$wildfish[i]
      Pars['F_w'] <- TR_TH$wildfish[i]
      
      M.0 = TR_TH$thresh_val[i]
      
      # redefine root & event with correct threshold and effecitivity
      rootfunc <- function(Time, y, Pars) {return (y[1] + y[2] - M.0*F_d.0)} # threshold hit
      eventfunc <- function(Time, y, Pars) {
        y[1] <- (1-ep_r.0)*y[1]
        y[2] <- (1-ep_s.0)*y[2]
        return(y)
      } # at the root, do this event
      
      # Pars['M.0']
      
      State = c(L_rd = 0, L_sd = 0.001, L_rj = 0, L_sj = 0, F_j = n_N1.0, F_a = n_N1.0, L_rn = 0, L_sn = n_N1.0, F_n = n_N1.0, L_rw = 0, L_sw = F_w.0)
      out <- ode(Time, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
      # then take output from out and use to start new simulation state
      State = out[dim(out)[1], 2:2:dim(out)[2]]; State["L_sd"] = State["L_sd"]-0.001; State["L_rd"] = 0.001 # can just replace a suseptable louse with a resistant louse
      Time_invade = seq(0, 5*n-31, length.out = 5*n-30)
      out <- ode(Time_invade, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
      Time_outs = seq(0, 30, length.out = 31) # simulate just last month
      State = out[dim(out)[1], 2:dim(out)[2]]
      out <- ode(Time_outs, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE, maxroots = 20000))
      
      out_df <- as.data.frame(out)
      # get the last 30 days
      # out_last_month <- out_df[(length(Time)-30):length(Time), ]
      
      resist_percent <- mean(out_df$L_rd/(out_df$L_rd + out_df$L_sd))
   
      }}

# end parallel cluster
stopCluster(cl)

# manipulate dataframes
{
evol_dat <- TR_TH
evol_dat$value <- unlist(resist_percent_list)
}

evol_outs <- list(base_pars, evol_dat)
# save(evol_outs, file = "simulated data/evol_outs.Rdata")

base_par_range <- evol_dat %>% filter(type == "Base parameters") %>% select(-type) %>% 
  pivot_wider(names_from = treat_eff, values_from = value) %>% 
  select(c(thresh_val, `0.3`, `0.9`))
colnames(base_par_range) <- c("thresh_val", "lower", "upper")

# make figures
{fig4 <- ggplot(data = evol_dat) + 
    geom_linerange(data = base_par_range,
                   aes(x = thresh_val, ymin = lower, ymax = upper), size = 1.4, col = "gray") + 
    geom_point(aes(thresh_val, value, col = as.factor(treat_eff)), size = 1.4) + 
    facet_wrap(~type, nrow=2) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
    theme(legend.position = "bottom",
          legend.title = element_text(size = 10), 
          legend.text = element_text(size = 10), 
          axis.title = element_text(size = 10), 
          axis.text = element_text(size = 10), 
          strip.text = element_text(size = 10)) +
    labs(x = "Treatment threshold", y = "Proportion resistant", col = "Treatment efficacy") 
}

fig4 # has some more efficacy options

#----figure 5: frontiers----

# need high treatment cost & link population size
{tmp <- merge(econ_dat[,c(1,2,6)], link_size[,c(4:6)], by = c("Treat_Int", "Threshold"))
colnames(tmp) <- c("treat_eff", "thresh_val", "econ", "ecol")
# need base parameter case
tmp2 <- evol_dat %>% filter(type == "Base parameters") %>% select(c("treat_eff", "thresh_val", "value"))
colnames(tmp2) <- c("treat_eff", "thresh_val", "evol")
fronts <- merge(tmp, tmp2, by = c("treat_eff", "thresh_val"))
}

# need to add lines
{ecol_econ <- ggplot(fronts, aes(ecol, -econ, col = as.factor(treat_eff))) + 
    # line across threshold options
    geom_path(aes(group = as.factor(thresh_val)), col = "gray80", lty = "dashed") + 
    # line across fronters
    geom_path(data = fronts[c(25:19, 10, 1), ], 
              col = "gray45") + 
    guides(lty = "none") + 
    geom_point(size = 1.4) + 
    scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10), 
        legend.position="none") +
  labs(x = "Juv. salmon pop. size", y = "Neg. economic costs", col = "Treatment \nefficacy") 
evol_econ <- ggplot(fronts, # %>% filter(thresh_val %in% M.vec[c(1, 3, 5, 7, 9, 11)]), 
                    aes(-evol, -econ)) + 
  # line across threshold options
  geom_path(aes(group = as.factor(thresh_val)), col = "gray80", lty = "dashed") + 
  # line across fronters 
  geom_line(data = fronts[c(21, 22, 23, 24, 15, 6, 8, 19, 10, 1), ], 
            col = "gray45") + 
  geom_point(size = 1.4, aes(col = as.factor(treat_eff))) + 
  scale_color_viridis(discrete = T, direction = -1) +  
  # scale_linetype_manual(values=c("1F", "dotted", "dotdash", "dashed")) + 
  theme_classic() + # adding threshold labels manually
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10), 
        legend.position="none") +
  labs(x = "Neg. prop. resist.", y = "Neg. economic costs", col = "Treatment \nefficacy") #, lty = "Treatment \nthreshold") 
# look at with threshold values??
ecol_evol <- ggplot(fronts, aes(ecol, -evol, col = as.factor(treat_eff))) + 
  # line across threshold options
  geom_path(aes(group = as.factor(thresh_val)), col = "gray80", lty = "dashed") + 
  # line across fronters
  geom_line(data = fronts[c(1:9, 16), ], 
            col = "gray45") + 
  geom_point(size = 1.4) + guides(lty = "none") + 
  scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10), 
        theme(legend.position="none")) +
  labs(x = "Juv. salmon pop. size", y = "Neg. prop. resist.", col = "Treatment \nefficacy") 
}

leg <- get_legend(ecol_evol + 
                    guides(color = guide_legend(nrow = 1)) +
                    theme(legend.position = "bottom", 
                          legend.box = "horizontal"))
plts_fig5 <- plot_grid(ecol_econ + theme(legend.position = "none"), 
                  ecol_evol + theme(legend.position = "none"),
                  evol_econ + theme(legend.position = "none"),
                  nrow = 1)

# fig5 <- plot_grid(plts_fig5, leg, 
#                   rel_heights = c(1, 1),
#                   nrow = 2)

fig5 <- plts_fig5

fig5

#----figure 6: trade-offs----

# need to smooth out and maybe increase the selectivity range??
# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

# set up & make function
{
cost.vect <- seq(from = 0.9, to = 0.99, length.out = 10)
select.vect <- seq(0.01, to = 0.1, length.out = 10) 
cost_select <- expand.grid(select = select.vect, cost = cost.vect)
  
  resist_cost_select <- foreach(i=1:dim(cost_select)[1], .packages="deSolve") %dopar% {
    
    ep_s.0 = 0.9 # using fixed set up for efficacy
    ep_r.0 = ep_s.0*cost_select$select[i] 
    l_r.0 = l_s.0*cost_select$cost[i] 
    
    Pars = base_pars
    Pars['ep_s'] <- ep_s.0
    Pars['ep_r'] <- ep_r.0
    Pars['l_r'] <- l_r.0
    
    M.0 = 0.7 # moderate values, in base case, produces resistance between 0.75-0.9ish
    
    # redefine root & event with correct threshold and effecitivity
    rootfunc <- function(Time, y, Pars) {return (y[1] + y[2] - M.0*F_d.0)} # threshold hit
    eventfunc <- function(Time, y, Pars) {
      y[1] <- (1-ep_r.0)*y[1]
      y[2] <- (1-ep_s.0)*y[2]
      return(y)
    } # at the root, do this event
    
    # Pars['M.0']
    
    State = c(L_rd = 0, L_sd = 0.001, L_rj = 0, L_sj = 0, F_j = n_N1.0, F_a = n_N1.0, L_rn = 0, L_sn = n_N1.0, F_n = n_N1.0, L_rw = 0, L_sw = F_w.0)
    out <- ode(Time, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
    # then take output from out and use to start new simulation state
    State = out[dim(out)[1], 2:2:dim(out)[2]]; State["L_sd"] = State["L_sd"]-0.001; State["L_rd"] = 0.001 # can just replace a suseptable louse with a resistant louse
    Time_invade = seq(0, 5*n-31, length.out = 5*n-30)
    out <- ode(Time_invade, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
    Time_outs = seq(0, 30, length.out = 31) # simulate just last month
    State = out[dim(out)[1], 2:dim(out)[2]]
    out <- ode(Time_outs, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE, maxroots = 20000))
    
    out_df <- as.data.frame(out)
    # get the last 30 days
    # out_last_month <- out_df[(length(Time)-30):length(Time), ]
    
    resist_percent <- mean(out_df$L_rd/(out_df$L_rd + out_df$L_sd))
    
  }}

# end parallel cluster
stopCluster(cl)

# manipulate dataframes
{
  cost_select$value <- unlist(resist_cost_select)
}

bio_trade_out <- list(base_pars, cost_select)
# save(bio_trade_out, file = "simulated data/bio_trade_out.Rdata")

# make figures
{fig6 <- ggplot(cost_select, aes(cost, select, col = value)) + 
    geom_point(size = 3) +    
    scale_color_viridis(direction = -1, limits = c(0.75, 1), breaks = c(0.75, 0.875, 1)) +  
    theme_classic() + 
    geom_point(data = data.frame(select = 0.05, cost = 0.95), col = "black", size = 5, pch = 1) + 
    theme(legend.position = "right",
          legend.title = element_text(size = 5), 
          legend.text = element_text(size = 5), 
          axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.key.size = unit(0.35, 'cm')) +
    labs(x = "Resistance cost", y = "Resistance benefit", col = "Prop. \nResistant") 
}

fig6 # has some more efficacy options

#----figure 7: real-ish parameters----

# should run base_parameters section 1 again...

# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

# set up & make function
{ep.vec <- seq(from = 0.3, to = 0.9, length.out = 3) # treatment effect --> might be interesting to also do with wild fish to see how threshold & wild fish interact
  M.vec <- seq(from = 0.01, to = 60.01, length.out = 11) # threshold
  TR_TH <- expand.grid(treat_eff = ep.vec, thresh_val = M.vec)

  # think about how to prep the foreach loop --> check how you did it in the evo rescue files?
  
  approx_pars_out <- foreach(i=1:dim(TR_TH)[1], .packages="deSolve", .combine="rbind") %dopar% {
    
    ep_s.0 = TR_TH$treat_eff[i]
    ep_r.0 = ep_s.0*0.5 # make 25% of ep_s
    
    Pars = approx_pars
    Pars['ep_s'] <- ep_s.0; Pars['ep_r'] <- ep_r.0

    M.0 = TR_TH$thresh_val[i]
    F_d.0 = Pars['F_d']
    F_w.0 = Pars['F_w']
    
    # redefine root & event with correct threshold and effecitivity
    rootfunc <- function(Time, y, Pars) {return (y[1] + y[2] - M.0*F_d.0)} # threshold hit
    eventfunc <- function(Time, y, Pars) {
      y[1] <- (1-ep_r.0)*y[1]
      y[2] <- (1-ep_s.0)*y[2]
      return(y)
    } # at the root, do this event
    
    State = c(L_rd = 0, L_sd = 0.001, L_rj = 0, L_sj = 0, F_j = n_N1.0*10000, F_a = n_N1.0*10000, L_rn = 0, L_sn = n_N1.0*10000, F_n = n_N1.0*10000, L_rw = 0, L_sw = F_w.0)
    out <- ode(Time, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
    # then take output from out and use to start new simulation state
    State = out[dim(out)[1], 2:2:dim(out)[2]]; State["L_sd"] = State["L_sd"]-0.001; State["L_rd"] = 0.001 # can just replace a suseptable louse with a resistant louse
    Time_invade = seq(0, 5*n-31, length.out = 5*n-30)
    out <- ode(Time_invade, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
    Time_outs = seq(0, 30, length.out = 31) # simulate just last month
    State = out[dim(out)[1], 2:dim(out)[2]]
    out <- ode(Time_outs, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE, maxroots = 20000))
    
    out_last_month <- as.data.frame(out)
    
    resist_percent <- mean(out_last_month$L_rd/(out_last_month$L_rd + out_last_month$L_sd))
    num_treats <- length(attributes(out)$troot)
    farm_burd <- mean((out_last_month$L_rd + out_last_month$L_sd)/F_d.0)
    link_pop <- mean(out_last_month$F_j)
    
    list(resist_percent, link_pop, num_treats, farm_burd)
    
  }}

# end parallel cluster
stopCluster(cl)

# manipulate dataframes
{
  approx_out <- cbind(TR_TH, unlist(approx_pars_out[,1]), 
                             unlist(approx_pars_out[,2]), 
                             unlist(approx_pars_out[,3]), 
                             unlist(approx_pars_out[,4]))
  colnames(approx_out) <- c(colnames(TR_TH), "Proportion resistant", "Link population size", "num_treat", "burd_farm")
  approx_out$`Economic losses` <- 2.5*approx_out$num_treat + 0.1*approx_out$burd_farm
  approx_out$`Econ. losses \n(standardized)` <- approx_out$`Economic losses`/max(approx_out$`Economic losses`)
  approx_out$`Juv. salmon pop. size \n(standardized)` <- approx_out$`Link population size`/max(approx_out$`Link population size`)
  approx_out_long <- pivot_longer(approx_out[, c(1:3, 8:9)], cols = 3:5, names_to = "output_type")
}

approx_outputs <- list(approx_pars, approx_out, approx_out_long)
# save(approx_outputs, file = "simulated data/approx_outputs.Rdata")

# make figures
{fig7 <- ggplot(data = approx_out_long) + 
    geom_point(aes(thresh_val, value, col = as.factor(treat_eff)), size = 1.4) + 
    facet_wrap(~output_type, nrow=1, scale = "free_y") +
    scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
    theme(legend.title = element_text(size = 10), 
          legend.text = element_text(size = 10), 
          axis.title = element_text(size = 10), 
          axis.text = element_text(size = 10), 
          strip.text = element_text(size = 10), 
          legend.position = "bottom") +
    labs(x = "Treatment threshold", y = NULL, col = "Treatment efficacy") 
}

fig7 # has some more efficacy options

#----figure 8: zoom on real-ish parameters----

# should run base_parameters section 1 again...

# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

# set up & make function
{ep.vec <- seq(from = 0.3, to = 0.9, length.out = 3) # treatment effect --> might be interesting to also do with wild fish to see how threshold & wild fish interact
  M.vec <- seq(from = 0.01, to = 4.01, length.out = 21) # threshold
  TR_TH <- expand.grid(treat_eff = ep.vec, thresh_val = M.vec)
  
  # think about how to prep the foreach loop --> check how you did it in the evo rescue files?
  
  approx_pars_out <- foreach(i=1:dim(TR_TH)[1], .packages="deSolve", .combine="rbind") %dopar% {
    
    ep_s.0 = TR_TH$treat_eff[i]
    ep_r.0 = ep_s.0*0.5 # make 25% of ep_s
    
    Pars = approx_pars
    Pars['ep_s'] <- ep_s.0; Pars['ep_r'] <- ep_r.0
    
    M.0 = TR_TH$thresh_val[i]
    F_d.0 = Pars['F_d']
    F_w.0 = Pars['F_w']
    
    # redefine root & event with correct threshold and effecitivity
    rootfunc <- function(Time, y, Pars) {return (y[1] + y[2] - M.0*F_d.0)} # threshold hit
    eventfunc <- function(Time, y, Pars) {
      y[1] <- (1-ep_r.0)*y[1]
      y[2] <- (1-ep_s.0)*y[2]
      return(y)
    } # at the root, do this event
    
    State = c(L_rd = 0, L_sd = 0.001, L_rj = 0, L_sj = 0, F_j = n_N1.0*10000, F_a = n_N1.0*10000, L_rn = 0, L_sn = n_N1.0*10000, F_s = n_N1.0*10000, L_rw = 0, L_sw = F_w.0)
    out <- ode(Time, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
    # then take output from out and use to start new simulation state
    State = out[dim(out)[1], 2:2:dim(out)[2]]; State["L_sd"] = State["L_sd"]-0.001; State["L_rd"] = 0.001 # can just replace a suseptable louse with a resistant louse
    Time_invade = seq(0, 5*n-31, length.out = 5*n-30)
    out <- ode(Time_invade, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))
    Time_outs = seq(0, 30, length.out = 31) # simulate just last month
    State = out[dim(out)[1], 2:dim(out)[2]]
    out <- ode(Time_outs, y = State, func = systemTREAT, parms = Pars, rootfun = rootfunc, events = list(func = eventfunc, root = TRUE, maxroots = 20000))
    
    out_last_month <- as.data.frame(out)
    
    resist_percent <- mean(out_last_month$L_rd/(out_last_month$L_rd + out_last_month$L_sd))
    num_treats <- length(attributes(out)$troot)
    farm_burd <- mean((out_last_month$L_rd + out_last_month$L_sd)/F_d.0)
    link_pop <- mean(out_last_month$F_j)
    
    list(resist_percent, link_pop, num_treats, farm_burd)
    
  }}

# end parallel cluster
stopCluster(cl)

# manipulate dataframes
{
  approx_out_zoom <- cbind(TR_TH, unlist(approx_pars_out[,1]), 
                      unlist(approx_pars_out[,2]), 
                      unlist(approx_pars_out[,3]), 
                      unlist(approx_pars_out[,4]))
  colnames(approx_out_zoom) <- c(colnames(TR_TH), "Proportion resistant", "Link population size", "num_treat", "burd_farm")
  approx_out_zoom$`Economic losses` <- 2.5*approx_out_zoom$num_treat + 0.1*approx_out_zoom$burd_farm
  approx_out_zoom$`Econ. losses \n(standardized)` <- approx_out_zoom$`Economic losses`/max(approx_out$`Economic losses`)
  approx_out_zoom$`Juv. salmon pop. size \n(standardized)` <- approx_out_zoom$`Link population size`/max(approx_out$`Link population size`)
  approx_out_zoom_long <- pivot_longer(approx_out_zoom[, c(1:3, 8:9)], cols = 3:5, names_to = "output_type")
}

approx_outputs_zoom <- list(approx_pars, approx_out_zoom, approx_out_zoom_long)
# save(approx_outputs_zoom, file = "simulated data/approx_outputs.Rdata")

approx_out_long$scale <- rep("out")
approx_out_zoom_long$scale <- rep("in")
approx_both <- rbind(approx_out_long, approx_out_zoom_long)

# make figures
{fig8a <- ggplot(data = approx_out_zoom_long) + 
    geom_point(aes(thresh_val, value, col = as.factor(treat_eff)), size = 1.4) + 
    facet_wrap(~output_type) + coord_cartesian(ylim = c(0, 1)) + 
    scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
    theme(legend.title = element_text(size = 10), 
          legend.text = element_text(size = 10), 
          axis.title = element_text(size = 10), 
          axis.text = element_text(size = 10), 
          strip.text = element_text(size = 10), 
          legend.position = "none") +
    labs(x = NULL, y = NULL, col = "Treatment efficacy") 
}

{fig8b <- ggplot(data = approx_out_long) + 
    geom_point(aes(thresh_val, value, col = as.factor(treat_eff)), size = 1.4) + 
    facet_wrap(~output_type) + coord_cartesian(ylim = c(0, 1)) +
    scale_color_viridis(discrete = T, direction = -1) +  theme_classic() + 
    theme(legend.title = element_text(size = 10), 
          legend.text = element_text(size = 10), 
          strip.text = element_blank(),
          axis.title = element_text(size = 10), 
          axis.text = element_text(size = 10), 
          legend.position = "bottom") +
    labs(x = "Treatment threshold", y = NULL, col = "Treatment efficacy") 
}

fig8 <- plot_grid(fig8a, fig8b, ncol = 1, rel_heights = c(1, 1.1))

fig8 # has some more efficacy options

#----print figures----

png("plots/economic_outcomes.png",height=170,width=85,res=400,units='mm')
print(fig2)
dev.off()

png("plots/ecological_outcomes.png",height=80,width=85,res=400,units='mm')
print(fig3)
dev.off()

png("plots/evolution_outcomes.png",height=180,width=170,res=400,units='mm')
print(fig4)
dev.off()

png("plots/expend_fronts.png",height=80,width=170,res=400,units='mm')
print(fig5)
dev.off()

png("plots/bio_tradeoff.png",height=80,width=85,res=400,units='mm')
print(fig6)
dev.off()

png("plots/approx_params_wide.png",height=80,width=170,res=400,units='mm')
print(fig7)
dev.off()

png("plots/approx_params.png",height=180,width=170,res=400,units='mm')
print(fig8)
dev.off()