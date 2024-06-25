# continuous dynamics
systemTREAT <- function (Time, y, Pars) {
    with(as.list(c(State, Pars)), {
       
        # domestic parasites
        dy1 = ifelse(y[1] <= 0, 0, 
                     (l_r/(v_r + B_d*F_d))*F_d*B_d*(p_dd*y[1] + p_dl*y[7]) - (u_r+Z)*y[1] - t_d*y[1]*(y[1] + y[2])/F_d)
        dy2 = ifelse(y[2] <= 0, 0, 
                     (l_s/(v_s + B_d*F_d))*F_d*B_d*(p_dd*y[2] + p_dl*y[8]) - (u_s+Z)*y[2] - t_d*y[2]*(y[1] + y[2])/F_d)
        
        # domesitc hosts in equilibirum
      
        # link juvenile parasites
        dy3 = (l_r/(v_r + B_l*y[5]))*y[5]*B_l*(p_ld*y[1]) - (u_r+z_j+T_out+alp)*y[3] - alp*y[3]*(y[3]+y[4])/y[5]
        dy4 = (l_s/(v_s + B_l*y[5]))*y[5]*B_l*(p_ld*y[2]) - (u_s+z_j+T_out+alp)*y[4] - alp*y[4]*(y[3]+y[4])/y[5]
        
        # link juvenile hosts
        dy5 = n_N1*y[9]/(n_N2+y[9]) - alp*(y[3] + y[4]) - z_j*y[5] - T_out*y[5] # hosts leave by maturing
        
        # link adult hosts
        dy6 = T_out*y[5] - T_in*y[6]
        
        # link spawning parasites
        dy7 = (l_r/(v_r + B_l*y[9]))*y[9]*B_l*p_ll*y[7] + T_in*y[10]*(y[6]/(y[6] + F_w)) - (u_r+z_n)*y[7] 
        dy8 = (l_s/(v_s + B_l*y[9]))*y[9]*B_l*p_ll*y[8] + T_in*y[11]*(y[6]/(y[6] + F_w)) - (u_r+z_n)*y[8] 
        
        # link spawning hosts
        dy9 = T_in*y[6] - z_n*y[9]
        
        # wild subsidy
        dy10 = (l_r/(v_r + B_w*(y[6] + F_w)))*B_w*(y[6] + F_w)*y[10] + T_out*y[3] - (z_w*(F_w/(F_w+y[6])) + u_r + T_in*(y[6]/(y[6] + F_w)))*y[10] - t_w*y[10]*(y[10]+y[11])/(F_w+y[6])
        dy11 = (l_s/(v_s + B_w*(y[6] + F_w)))*B_w*(y[6] + F_w)*y[11] + T_out*y[4] - (z_w*(F_w/(F_w+y[6])) + u_s + T_in*(y[6]/(y[6] + F_w)))*y[11] - t_w*y[11]*(y[10]+y[11])/(F_w+y[6])
        
        # wild hosts in equilibrium
        
        derivs = c(dy1, dy2, dy3, dy4, dy5, dy6, dy7, dy8, dy9, dy10, dy11)
        
        derivs[!is.finite(derivs)] = 0
        
        return(list(derivs))
    })
}


# treat at discrete time points
rootfunc <- function(Time, y, Pars) {return (y[1] + y[2] - th.0*F_d.0)}

# how effective is treatment?
eventfunc <- function(Time, y, Pars) {
  y[1] <- (1-t_r.0)*y[1]
  y[2] <- (1-t_s.0)*y[2]
  return(y)
}

