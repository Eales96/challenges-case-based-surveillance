

SIRS <- function(t, beta, gamma, delta,S0, I0, R0){
  out <- data.frame(time=t,
                    S=S0,
                    I=I0,
                    R=R0,
                    inc=0)
  
  for(i in 2:nrow(out)){
    S = out$S[i-1]
    I = out$I[i-1]
    R = out$R[i-1]
    
    Sd <- delta*R  - beta*S*I
    Id <- beta*S*I - gamma*I
    Rd <- gamma*I  - delta*R
    
    out$S[i] <- S +Sd
    out$I[i] <- I +Id
    out$R[i] <- R +Rd
    out$inc[i] <- beta*S*I
    
  }
  
  out
}



calculate_time_series <- function(I_t, S_t, Ts_t, Ta_t, F_tau, P_tau, epsilon){
  
  out <- data.frame(time = seq(1,length(I_t),1),
                    I_t = I_t,
                    S_t = S_t,
                    Ts_t = Ts_t,
                    Ta_t = Ta_t,
                    ns_t = NA,
                    cs_t1 = NA,
                    cs_t2 = NA,
                    cs_t3 = NA,
                    ca_t1 = NA,
                    ca_t2 = NA,
                    cs_t = NA,
                    ca_t = NA)
  
  for(i in 31:nrow(out)){
    
    
    integral0 <- sum(rev(out$I_t[(i-30):i]) * F_tau)
    integral1 <- sum(rev(out$I_t[(i-30):i]) * F_tau * P_tau)
    integral2 <- sum(rev(out$I_t[(i-30):i]) * (1-F_tau) * P_tau)
    
    out$ns_t[i] <- integral0 + out$S_t[i]*(1 - integral0)
    out$cs_t1[i] <- out$Ts_t[i] * integral1
    out$cs_t2[i] <- out$Ts_t[i] * out$S_t[i] * integral2
    out$cs_t3[i] <- out$Ts_t[i] * epsilon * out$ns_t[i]
    
    out$ca_t1[i] <- out$Ta_t[i] * (1-out$S_t[i]) * integral2
    out$ca_t2[i] <- out$Ta_t[i] * epsilon * (1-out$ns_t[i])
    
    out$cs_t[i] <- out$cs_t1[i] + out$cs_t2[i] + out$cs_t3[i]
    out$ca_t[i] <- out$ca_t1[i] + out$ca_t2[i]
  }
  
  out
  
}