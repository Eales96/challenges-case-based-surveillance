#setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/BiasesInCaseData")

source('R/functions.R')
library(ggplot2)
library(cowplot)
##################################################################################
# Functions 
P_tau_calc <- function(tau, beta1, beta2, beta3, C){
  x = tau - C
  val = beta1 + x*beta2 +x*beta2*beta3*(x>0)
  exp(val)/(1+exp(val))
}

wt_sum <- function(pdf1, pdf2){
  dist <- rep(0, length(pdf1)+length(pdf2)-1)
  for(i in seq_len(length(pdf1))){
    dist[i:(i-1+length(pdf2))] <- dist[i:(i-1+length(pdf2))] + pdf1[i]*pdf2
  }
  dist*0.8/sum(pdf1)
}

F_tau_calc <- function(tau, shape1, rate1, shape2, rate2){
  
  pdf1 <- 0.8 * dgamma(tau, shape=shape1, rate=rate1)
  cdf2 <- pgamma(tau, shape=shape2, rate=rate2)
  
  
  ftau <- wt_sum(pdf1,
                 1-cdf2)
  
  ftau
  
}
########################################################################################


t <- seq(0,400,1)
test <- SIRS(t, beta=2., gamma=0.5, delta=0.01, S0=0.5-0.001, I0=0.001, R0=0.5 )
plot(test$time, test$inc)

I_t <- test$inc



S_t1 <- rep(0.0, nrow(test))
S_t2 <- rep(0.1, nrow(test))
t<- seq(0,350,1)
S_t3 <-c(rep(0,50), 0.7*(0.1-0.05*cos(2*pi*t/200)+0.0005*t-0.05*cos(2*pi*t/400)+0.05*sin(2*pi*t/100)))

Ts_t <- rep(1, nrow(test))
Ta_t <- rep(1, nrow(test))

shp1 <- 6**2 / 3.1**2
rt1 <- 6/ (3.1**2)

shp2 <- 8**2/4.84**2
rt2 <- 8/(4.84**2)

f_tau<-F_tau_calc(seq(0,15,1), shape1=9, rate1=3, shape2= 40, rate2=5)
#f_tau<-F_tau_calc(seq(0,15,1), shape1=shp1, rate1=rt1, shape2= shp2, rate2=rt2)
p_tau <- P_tau_calc(seq(0,30,1), beta1=1.51, beta2=2.19, beta3 = -1.1, C=3.18)


num <- sum(f_tau * p_tau)
den <- sum((1-f_tau) * p_tau)

num/den

epsilon <- 0


ts <- data.frame()
for(i in seq_len(10)){
  S_t <- rep((i-1)*0.02, nrow(test))
  print(i)
  
  ts_tmp<-calculate_time_series(I_t, S_t, Ts_t, Ta_t,
                                f_tau, p_tau, epsilon)
  ts_tmp$tag <- (i-1)*0.02
  ts <- rbind(ts, ts_tmp)
  
}

ts <-ts[ts$time>=50 & ts$time<250,]


###########################################################################################
###########################################################################################

scale1 <- max(ts$I_t)/max(ts$cs_t/ts$ns_t)
scale2 <- max(ts$I_t)/max(ts$ca_t/(1-ts$ns_t))
scale3 <- max(ts$I_t)/max(ts$ca_t+ts$cs_t)
scale4 <- max(ts$I_t)/max((5*ts$ca_t+ts$cs_t)/(5-4*ts$ns_t))
scale5 <- max(ts$I_t)/max((ts$ca_t+5*ts$cs_t)/(1+4*ts$ns_t))
scale6 <- max(ts$I_t)/max((ts$ca_t+10*ts$cs_t)/(1+9*ts$ns_t))


scale1b <- max(ts$cs_t/ts$ns_t)
scale2b <- max(ts$ca_t/(1-ts$ns_t))
scale3b <- max(ts$ca_t+ts$cs_t)
scale4b <- max((5*ts$ca_t+ts$cs_t)/(5-4*ts$ns_t))
scale5b <- max((ts$ca_t+5*ts$cs_t)/(1+4*ts$ns_t))
scale6b <- max((ts$ca_t+10*ts$cs_t)/(1+9*ts$ns_t))


plt1<-ggplot(ts)+
  geom_line(aes(x=time, y=cs_t/ns_t, color=factor(tag)))+
  geom_line(aes(x=time, y=I_t/scale1), linetype='dashed', color='red3')+
  scale_color_viridis_d("S(t)=", option='viridis')+
  scale_y_continuous(expression("Proportion of tests positive"), sec.axis = sec_axis(~.*scale1, "Incidence, I(t)"))+
  xlab("Time")+
  theme_bw()+
  geom_label(aes(x = 150, y = scale1b*1.05, label = "Symptomatic testing"), fill = "white")+
  theme(legend.position="none")


plt2<-ggplot(ts[ts$S_t==0,])+
  geom_line(aes(x=time, y=ca_t/(1-ns_t), color=factor(tag)))+
  geom_line(aes(x=time, y=I_t/scale2), linetype='dashed', color='red3')+
  scale_color_viridis_d("S(t)=", option='viridis')+
  scale_y_continuous(expression("Proportion of tests positive"), sec.axis = sec_axis(~.*scale2, "Incidence, I(t)"))+
  xlab("Time")+
  theme_bw()+
  geom_label(aes(x = 150, y = scale2b*1.05, label = "Asymptomatic testing"), fill = "white")+
  theme(legend.position="none",
        legend.background = element_rect(color='black'))


plt2leg<-ggplot(ts)+
  geom_line(aes(x=time, y=ca_t/(1-ns_t), color=factor(tag)))+
  geom_line(aes(x=time, y=I_t/scale2), linetype='dashed', color='red3')+
  scale_color_viridis_d("S(t)=", option='viridis')+
  scale_y_continuous(expression("Proportion of tests positive"), sec.axis = sec_axis(~.*scale2, "Incidence, I(t)"))+
  xlab("Time")+
  theme_bw()+
  geom_label(aes(x = 150, y = scale2b*1.05, label = "Asymptomatic testing"), fill = "white")+
  theme(legend.position="bottom",
        legend.background = element_rect(color='black'))+
  guides(color=guide_legend(nrow=1))


plt3<-ggplot(ts[ts$S_t==0,])+
  geom_line(aes(x=time, y=(cs_t+ca_t)/(1), color=factor(tag)))+
  geom_line(aes(x=time, y=I_t/scale3), linetype='dashed', color='red3')+
  scale_color_viridis_d("S(t)=", option='viridis')+
  scale_y_continuous(expression("Proportion of tests positive "), sec.axis = sec_axis(~.*scale3, "Incidence, I(t)"))+
  xlab("Time")+
  theme_bw()+
  geom_label(aes(x = 150, y = scale3b*1.05, label = "Random sampling"), fill = "white")+
  theme(legend.position="none",
        legend.background = element_rect(color='black'))

plt4<-ggplot(ts)+
  geom_line(aes(x=time, y=(cs_t+5*ca_t)/(5-4*ns_t), color=factor(tag)))+
  geom_line(aes(x=time, y=I_t/scale4), linetype='dashed', color='red3')+
  scale_color_viridis_d("S(t)=", option='viridis')+
  scale_y_continuous(expression("Proportion of tests positive "), sec.axis = sec_axis(~.*scale4, "Incidence, I(t)"))+
  xlab("Time")+
  theme_bw()+
  geom_label(aes(x = 150, y = scale4b*1.05, label = "T[a](t)==5*T[s](t)"),parse = T, fill = "white")+
  theme(legend.position="none",
        legend.background = element_rect(color='black'))

plt5<-ggplot(ts)+
  geom_line(aes(x=time, y=(5*cs_t+ca_t)/(4*ns_t+1), color=factor(tag)))+
  geom_line(aes(x=time, y=I_t/scale5), linetype='dashed', color='red3')+
  scale_color_viridis_d("S(t)=", option='viridis')+
  scale_y_continuous(expression("Proportion of tests positive "), sec.axis = sec_axis(~.*scale5, "Incidence, I(t)"))+
  xlab("Time")+
  theme_bw()+
  geom_label(aes(x = 150, y = scale5b*1.05, label = "T[a](t)==0.2*T[s](t)"),parse=T, fill = "white")+
  theme(legend.position="none",
        legend.background = element_rect(color='black'))

plt6<-ggplot(ts)+
  geom_line(aes(x=time, y=(10*cs_t+ca_t)/(9*ns_t+1), color=factor(tag)))+
  geom_line(aes(x=time, y=I_t/scale6), linetype='dashed', color='red3')+
  scale_color_viridis_d("S(t)=", option='viridis')+
  scale_y_continuous(expression("Proportion of tests positive "), sec.axis = sec_axis(~.*scale6, "Incidence, I(t)"))+
  xlab("Time")+
  theme_bw()+
  geom_label(aes(x = 150, y = scale6b*1.05, label = "T[a](t)==0.1*T[s](t)"),parse=T, fill = "white")+
  theme(legend.position="none",
        legend.background = element_rect(color='black'))



plt1<- plt1+theme(axis.title = element_blank())
plt2<- plt2+theme(axis.title = element_blank())
plt3<- plt3+theme(axis.title = element_blank())
plt4<- plt4+theme(axis.title = element_blank())
plt5<- plt5+theme(axis.title = element_blank())
plt6<- plt6+theme(axis.title = element_blank())


final_grid <- plot_grid(plt1,plt6, plt5, plt3, plt4, plt2, nrow=2)


library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(cowplot)
library(ggpubr)
leg <- get_legend(plt2leg)


y.grob <- textGrob("Proportion of tests positive", 
                   gp=gpar(fontsize=11), rot=90)

y.grob2 <- textGrob("Incidence", 
                   gp=gpar(fontsize=11), rot=270)

x.grob <- textGrob("Time", 
                   gp=gpar(fontsize=11))

sim_grd1<-arrangeGrob(final_grid, left = y.grob, bottom = x.grob, right=y.grob2)
sim_grd <- arrangeGrob(sim_grd1, bottom=leg)
sim_grd<-ggdraw()+draw_grob(sim_grd)


ggsave('figures/figure6.pdf', width=12, height=8, units= 'in')

################################################################################################
# Supplementary
###########################################################################################
