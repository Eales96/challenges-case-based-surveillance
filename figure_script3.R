#setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/BiasesInCaseData")

source('R/functions.R')
library(ggplot2)




# Function to return F(tau)*P(tau)

calc_ratio <- function(tf, hf, tp, hp){
  tau <- seq(1, 30, 1)
  f_tau <- rep(hf, 30) * (tau<=tf)
  p_tau <- rep(hp, 30) * (tau<=tp)
  
  num <- sum(f_tau * p_tau)
  den <- sum((1-f_tau) * p_tau)
  
  return(num/den)
}



tf_list <- seq(1, 20, 1)
tp_list <- seq(1, 20, 1)
hf_list <- seq(0.05, 0.95, length=19)
hp_list <- seq(0.5, 1, length=1)

df<- expand.grid(tf_list, tp_list,
                 hf_list, hp_list)

colnames(df) <- c("tf", "tp", "hf", "hp")

df$ratio <- -99
for(i in seq_len(nrow(df))){
  df$ratio[i] <- calc_ratio(df$tf[i], df$hf[i], df$tp[i], df$hp[i])
}


df$ratio[10000:20000]


plt1m<-ggplot(df, aes(x=tp/tf, y=ratio, color=factor(hf)))+
  geom_line()+
  #coord_cartesian(ylim=c(0,10))
  scale_y_log10()+
  coord_cartesian(xlim=c(0,20))+
  scale_color_viridis_d("Proportion symptomatic",option='turbo')+
  theme_bw()+
  xlab("Duration test positive / duration exhibiting symptoms")+
  ylab(label =expression("Ratio of symptomatic cases to asymptomatic cases (factor of "~T[s]/T[a]~ ")" ))+
  annotation_logticks(base=10, sides='l')+
  geom_vline(xintercept = 1, color='black', linetype="dashed")+
  guides(color=guide_legend(ncol=3,byrow=TRUE))+
  theme(legend.position = c(0.8,0.75),
        legend.background=element_rect(color='black'),
        panel.grid.minor.y = element_blank())



plt1i<-ggplot(df[df$tp==1 &df$tf==2,], aes(x=hf, y=ratio, color=factor(hf)))+
  geom_line(aes(x=hf, y=ratio), color='black', linetype="dashed")+
  geom_point(size=1)+
  scale_y_log10()+
  scale_x_continuous(breaks=c(0.0,0.2, 0.4, 0.6, 0.8,1.0))+
  scale_color_viridis_d(option='turbo')+
  theme_bw(base_size = 8)+
  xlab("Proportion symptomatic")+
  ylab(label =expression(~"Ratio of symptomatic cases\nto asymptomatic cases" ))+
  annotation_logticks(base=10, sides='l')+
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank(),
        plot.background = element_rect(color='black'),
        plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm"))


library(cowplot)

plt1 <- ggdraw()+
  draw_plot(plt1m)+
  draw_plot(plt1i, x=0.25, y=0.62, width=0.35, height=0.35)+
  labs(tag = "B")+
  theme(plot.tag.position = c(0.02,0.98))
plt1

ggsave('figures/figure3.pdf', width=7, height=6, units= 'in')


###################################################################################################

t <- seq(0,400,1)
test <- SIRS(t, beta=1, gamma=0.5, delta=0.015, S0=0.999, I0=0.000001, R0=0 )
plot(test$time, test$inc)

I_t <- test$inc
S_t <- rep(0, nrow(test))
Ts_t <- rep(1, nrow(test))
Ta_t <- rep(1, nrow(test))


I_t <- c(rep(0.0000001,50), test$inc)
S_t <- rep(0, length(I_t))
Ts_t <- rep(1, length(I_t))
Ta_t <- rep(1, length(I_t))

F_tau1 <- rep(0.8, 31) * (seq(0,30,1)<=5)
F_tau2 <- rep(0.75, 31) * (seq(0,30,1)<=8)
P_tau <- rep(1, 31) * (seq(0,30,1)<=10)


epsilon <- 0



#F_tau1 <- F_tau1[1:11]
#F_tau2 <- F_tau2[1:11]
#P_tau <- P_tau[1:11]


ts1<-calculate_time_series(I_t, S_t, Ts_t, Ta_t,
                           F_tau1, P_tau, epsilon)

ts2<-calculate_time_series(I_t, S_t, Ts_t, Ta_t,
                           F_tau2, P_tau, epsilon)



ts1$cs_t[1:61]<-NA
ts1$ca_t[1:61]<-NA

ts2$cs_t[1:61]<-NA
ts2$ca_t[1:61]<-NA

ts1$time <- ts1$time -50
ts2$time <- ts2$time -50
ts1 <- ts1[ts1$time >= 0,]
ts2 <- ts2[ts2$time >= 0,]

plt2m<-ggplot(ts1)+
  geom_line(aes(x=time, y=I_t))+
  geom_line(aes(x=time, y=0.02*cs_t/ca_t), color='red2')+
  geom_line(data=ts2, aes(x=time, y=0.02*cs_t/ca_t), color='blue2')+
  scale_y_continuous("Daily infection incidence", sec.axis = sec_axis(~./0.02, expression("Ratio of symptomatic cases to asymptomatic cases (factor of "~T[s]/T[a]~ ")" )))+
  xlab("Time")+
  theme_bw()+
  theme(legend.position="none",
        legend.backgroun=element_rect(color='black'))




F_tau1 <- rep(0.8, 3001) * (seq(0,30,0.01)<=5)
F_tau2 <- rep(0.75, 3001) * (seq(0,30,0.01)<=8)
P_tau <- rep(1, 3001) * (seq(0,30,0.01)<=10)

conv1 <- data.frame(tau = seq(0,30, 0.01),
                    F_tau1 = F_tau1,
                    F_tau2 = F_tau2,
                    P_tau = P_tau)

plt2i<-ggplot(conv1)+
  geom_line(aes(x=tau, y=F_tau1, color='F1(tau)'))+
  geom_line(aes(x=tau, y=F_tau2, color='F2(tau)'))+
  geom_line(aes(x=tau, y=P_tau, color='P(tau)'))+
  xlab(expression("Time since infection ("~tau~")"))+
  ylab("Probability")+
  theme_bw()+
  coord_cartesian(xlim=c(0,20))+
  scale_color_manual("Distribution",values = c('red2', 'blue2', 'green4'), labels = c(expression(F[1](tau)), expression(F[2](tau)),expression(P(tau))))+
  theme(legend.position=c(0.77,0.6),
        legend.backgroun=element_rect(color='black'),
        plot.background = element_rect(color='black'))


leg_tit <- expression("Value of "~frac(integral(I(t-tau)*F(tau)*P(tau)*d~tau, 0, infinity), integral(I(t-tau)*(1-F(tau))*P(tau)*d~tau, 0, infinity))~" for:")
plt2m<-ggplot(ts1)+
  geom_line(aes(x=time, y=I_t))+
  geom_line(aes(x=time, y=0.02*cs_t/ca_t, color='1'))+
  geom_line(data=ts2, aes(x=time, y=0.02*cs_t/ca_t, color="2"))+
  scale_y_continuous("Daily infection incidence", sec.axis = sec_axis(~./0.02, expression("Ratio of symptomatic cases to asymptomatic cases (factor of "~T[s]/T[a]~ ")" )))+
  xlab("Time")+
  theme_bw()+
  scale_color_manual(leg_tit,values = c('red2', 'blue2'), labels = c(expression(F[1](tau)), expression(F[2](tau))))+
  theme(legend.position=c(0.687,0.92),
        legend.background=element_rect(color='black'),
        legend.direction = 'horizontal')+
  guides(color=guide_legend(nrow=2,byrow=TRUE))


plt2 <- ggdraw()+
  draw_plot(plt2m)+
  draw_plot(plt2i, x=0.52, y=0.5, width=0.4, height=0.35)+
  labs(tag = "A")+
  theme(plot.tag.position = c(0.02,0.98))
plt2




plot_grid(plt2, plt1, nrow=1)
ggsave('figures/figure3.pdf', width=14, height=6, units= 'in')


