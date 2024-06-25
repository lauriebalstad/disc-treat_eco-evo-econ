#----system-specific pars----
# lice bio
l_s.0 = 6.35; l_r.0 = l_s.0*0.95 # lice reproduction -- note cost to resistence
v_s.0 = 1/5; v_r.0 = v_s.0 # juv. lice mortality
u_s.0 = 1/60; u_r.0 = u_s.0 # adult lice mortality
t_d.0 = 0.001; t_w.0 = 0.01 # density dep. on adult fish

# connectivity & attachement
p_dd.0 = 0.9 # stay on farm
p_ll.0 = p_dd.0 # stay in link
p_dl.0 = 1-p_dd.0 # spillover
p_ld.0 = p_dl.0 # spillback
p_ww.0 = 1 # all lice in wild stay in that environment, movement to link environment via migration
B_d.0 = 5.9*10^-10/1.3 # farm lice attachement
B_l.0 = B_d.0*0.75 # link lice attachement slightly lower in the link population -- bateman vibe
B_w.0 = B_d.0*0.25 # wild lice attachement -- bateman vibe

# treatment
ep_s.0 = 0.9 # suseptible lice effecivity
ep_r.0 = ep_s.0*0.5 # resistant lice effectivity

# host things
Z.0 = 0.67/365 # harvest
F_d.0 = 6*10^6 # normalize farm fish population size --> might need to be closer to 2e6?? for about a single farm v. whole network
F_w.0 = 20*F_d.0
n_N1.0 = 3*10^2
n_N2.0 = 2*10^2 
alp.0 = 0.02
z_w.0 = 0.005 # between harvest and juvenile salmon death
T_out.0 = 1/(0.25*365)
z_j.0 = T_out.0
T_in.0 = 1/(1.25*365)
z_n.0 = 1/(0.25*365)

n = 365*5 # for daily pars --> 10 years
Time = seq(0, n, length.out = n+1)

approx_pars = c(l_s = l_s.0, v_s = v_s.0, l_r = l_r.0, v_r = v_r.0, 
              u_s = u_s.0, u_r = u_r.0, t_d = t_d.0, t_w = t_w.0,
              p_dd = p_dd.0, p_ll = p_ll.0, p_dl = p_dl.0, p_ld = p_ld.0, p_ww = p_ww.0, 
              B_d = B_d.0, B_l = B_l.0, B_w = B_w.0,
              ep_r = ep_r.0, ep_s = ep_s.0, Z = Z.0,
              F_d = F_d.0, F_w = F_w.0, 
              n_N1 = n_N1.0, n_N2 = n_N2.0, alp = alp.0, z_j = z_j.0, z_w = z_w.0, z_n = z_n.0,
              T_out = T_out.0, T_in = T_in.0)

#----generic pars----
l_s.0 = 1; l_r.0 = 0.95*l_s.0
v_s.0 = 1 ; v_r.0 = v_s.0 
u_s.0 = .1; u_r.0 = .1 
t_d.0 = 0.15; t_w.0 = 0.15

p_dd.0 = 0.9 # low connectivity case (p_dd.0 close to 1)
p_ll.0 = p_dd.0 # lice in farm stay on farm...
p_dl.0 = 1-p_dd.0 # ...or jump to wild (spillover)
p_ld.0 = p_dl.0 # spillback rate
p_ww.0 = 1 # all lice in wild stay in that environment
B_d.0 = 1 # all the lice attach type vibe
B_l.0 = 0.75 # slightly lower in the link population -- bateman vibe
B_w.0 = 0.25 # lowest in wild

ep_s.0 = 0.9
ep_r.0 = ep_s.0*0.05 # similar to Murray

Z.0 = 0.2
F_d.0 = 1 # normalize farm fish population size
F_w.0 = 5 # I think this is right?
n_N1.0 = 1
n_N2.0 = 1
alp.0 = 0.4 # parasite induced mortality
z_j.0 = 0.1 # 0.1
z_w.0 = 0.1
T_out.0 = 0.2
T_in.0 = 0.1
z_n.0 = 0.5

n = 365 
Time = seq(0, 5*n, length.out = 5*n+1)

base_pars = c(l_s = l_s.0, v_s = v_s.0, l_r = l_r.0, v_r = v_r.0, 
         u_s = u_s.0, u_r = u_r.0, t_d = t_d.0, t_w = t_w.0,
         p_dd = p_dd.0, p_ll = p_ll.0, p_dl = p_dl.0, p_ld = p_ld.0, p_ww = p_ww.0, 
         B_d = B_d.0, B_l = B_l.0, B_w = B_w.0,
         ep_r = ep_r.0, ep_s = ep_s.0, Z = Z.0,
         F_d = F_d.0, F_w = F_w.0, 
         n_N1 = n_N1.0, n_N2 = n_N2.0, alp = alp.0, z_j = z_j.0, z_w = z_w.0, z_n = z_n.0,
         T_out = T_out.0, T_in = T_in.0)