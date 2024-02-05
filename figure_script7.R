setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/BiasesInCaseData")

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
  dist*0.4/sum(pdf1)
} # 0.8

F_tau_calc <- function(tau, shape1, rate1, shape2, rate2){
  
  pdf1 <- 0.8 * dgamma(tau, shape=shape1, rate=rate1)
  cdf2 <- pgamma(tau, shape=shape2, rate=rate2)
  
  
  ftau <- wt_sum(pdf1,
                 1-cdf2)
  
  ftau
  
}
########################################################################################


t <- seq(0,200,1)
#test1 <- SIRS(t, beta=0.8, gamma=0.16, delta=0.00005, S0=0.5-0.001, I0=0.0001, R0=0.5 )
#test <- SIRS(t, beta=0.8, gamma=0.16, delta=0.008, S0=0.5-0.001, I0=0.001, R0=0.5 )

test1 <- SIRS(t, beta=0.9, gamma=0.16, delta=0.00005, S0=0.5-0.0001, I0=0.0001, R0=0.5 )
test <- SIRS(t, beta=0.7, gamma=0.16, delta=0.005, S0=0.5-0.001, I0=0.001, R0=0.5 )

plot(test$time, test$inc)
plot(test1$time, test1$inc)

I_t <- c(rep(0, 49),test$inc)



S_t1 <- c(rep(0, 49),test1$inc*4)+0.005 #5
S_t1_alt <- c(rep(0, 49),test1$inc*0)+0.00 #5

Ts_t <- rep(1, 250)
Ta_t <- rep(1, 250)

f_tau<-F_tau_calc(seq(0,15,1), shape1=9, rate1=3, shape2= 40, rate2=10) #rate 2 =
p_tau <- P_tau_calc(seq(0,30,1), beta1=1.51, beta2=2.19, beta3 = -1.1, C=3.18)


num <- sum(f_tau * p_tau)
den <- sum((1-f_tau) * p_tau)

num/den

epsilon <- 0

ts<-calculate_time_series(I_t, S_t1, Ts_t, Ta_t,
                          f_tau, p_tau, epsilon)

ts_alt<-calculate_time_series(I_t, S_t1_alt, Ts_t, Ta_t,
                          f_tau, p_tau, epsilon)


ts <-ts[ts$time>50 & ts$time<250,]
ts_alt <-ts_alt[ts_alt$time>50 & ts_alt$time<250,]

plot(ts$time, ts$I_t)
plot(ts$time, ts$S_t)

plot(ts$time, ts$S_t)
plot(ts$time, ts$cs_t/ts$ns_t)
plot(ts$time, ts$ca_t/(1-ts$ns_t))
plot(ts$time, (ts$ca_t+ts$cs_t))

ts1<-ts
ts1$tag <- "Symptomatic testing"
ts1$prop <- ts1$cs_t/ts1$ns_t

ts2<-ts
ts2$tag <- "Asymptomatic testing"
ts2$prop <- ts2$ca_t/(1-ts2$ns_t)

ts3<-ts
ts3$tag <- "Random testing"
ts3$prop <- ts3$cs_t + ts3$ca_t
df <- rbind(ts1, ts2, ts3)

ggplot(data = ts)+
  geom_line(aes(x=time, y=I_t,color="Infection incidence"))+
  geom_line(aes(x=time, y=S_t,color="Background proportion symptomatic"))+
  geom_line(aes(x=time, y=cs_t/ns_t,color="Symptomatic testing"))+
  geom_line(aes(x=time, y=cs_t/(1-ns_t),color="Asymptomatic testing"))+
  geom_line(aes(x=time, y=(cs_t+ca_t),color="Random testing"))+
  theme_bw()


scale1 <- max(df[df$tag=="Symptomatic testing",]$I_t) /max( df[df$tag=="Symptomatic testing",]$prop)
scale2 <- max(df[df$tag=="Asymptomatic testing",]$I_t) /max( df[df$tag=="Asymptomatic testing",]$prop)
scale3 <- max(df[df$tag=="Random testing",]$I_t) /max( df[df$tag=="Random testing",]$prop)

df_new <- df

df$fudge <- rep(ts_alt$ns_t,3)
df[df$tag=="Symptomatic testing",]$tag <- "ILI testing"
df[df$tag=="Random testing",]$tag <- "Random population testing"



# influenza - fudge
# sars-cov-2 (S_t-0.005)*(1-fudge)
#(0.005)*(1-fudge)

plt4<-ggplot(data = df[df$tag=="ILI testing",])+
  geom_area(aes(x=time-50, y=S_t*(1-fudge)+fudge, fill="SARS-CoV-2"))+
  geom_area(aes(x=time-50, y=0.005*(1-fudge) + fudge, fill="Influenza"))+
  geom_area(aes(x=time-50, y=0.005*(1-fudge), fill="Allergies"))+
  geom_line(aes(x=time-50, y=ns_t))+
  scale_fill_brewer("Contribution to ILI", palette="Dark2")+
  geom_line(aes(x=time-50, y=I_t*6.1),color="black", linetype="dashed")+
  xlab("Time")+
  scale_y_continuous("Proportion of population currently exhibiting ILI", sec.axis = sec_axis(~.*1/6.1, "Infection incidence"))+
  theme_bw()+
  theme(legend.position = c(0.78,0.78),
        legend.background = element_rect(color='black'))



plt1<-ggplot(data = df[df$tag=="ILI testing",])+
  geom_line(aes(x=time-50, y=prop, color="Proportion positive\n(observed)"))+
  geom_line(aes(x=time-50, y=I_t/scale1, color="Infection incidence\n(unobserved)"), linetype="dashed")+
  facet_wrap(.~tag)+
  xlab("Time")+
  scale_color_manual(values=c('black','red2'))+
  scale_y_continuous("Proportion of tests positive for influenza", sec.axis = sec_axis(~.*scale1, "Infection incidence"))+
  theme_bw()

plt2<-ggplot(data = df[df$tag=="Asymptomatic testing",])+
  geom_line(aes(x=time-50, y=prop, color="Proportion positive\n(observed)", linetype="Proportion positive"))+
  geom_line(aes(x=time-50, y=I_t/scale2, color="Infection incidence\n(unobserved)", linetype="Infection incidence"))+
  facet_wrap(.~tag)+
  xlab("Time")+
  scale_color_manual("",values=c('black','red2'))+
  scale_linetype_manual("",values=c("dashed","solid"))+
  scale_y_continuous("Proportion of tests positive", sec.axis = sec_axis(~.*scale2, "Infection incidence"))+
  theme_bw()

plt3<-ggplot(data = df[df$tag=="Random population testing",])+
  geom_line(aes(x=time-50, y=prop, color="Proportion positive\n(observed)",linetype="Proportion positive\n(observed)"))+
  geom_line(aes(x=time-50, y=I_t/scale3, color="Infection incidence\n(unobserved)", linetype="Infection incidence\n(unobserved)"))+
  facet_wrap(.~tag)+
  xlab("Time")+
  scale_color_manual(name="test",values=c('black','red2'))+
  scale_linetype_manual(name="test",values=c('dashed','solid'))+
  scale_y_continuous("Proportion of tests positive for influenza (infection prevalence)", sec.axis = sec_axis(~.*scale3, "Infection incidence"))+
  theme_bw()


plt1<-plt1 + theme(legend.position = "none")+labs(tag="B")+
  theme(plot.tag.position = c(0.02,0.98))

plt3<-plt3 + theme(legend.position = c(0.75,0.75),
                   legend.background = element_rect(color = 'black'),
                   legend.title = element_blank())+labs(tag="C")+
  theme(plot.tag.position = c(0.02,0.98))

plt4<-plt4+labs(tag="A")+
  theme(plot.tag.position = c(0.02,0.98))


library(grid)
library(gridExtra)
final_grid <- plot_grid(plt4, plt1, plt3, nrow=1)
ggsave('figures/figure7.pdf', width=14, height=6, units= 'in')




ggplot(data = df[df$tag=="ILI testing",])+
  geom_line(aes(x=time-50, y=prop, color="Proportion positive\n(observed)"))+
  geom_line(aes(x=time-50, y=prop*cs_t*30, color="Proportion positive\n(observed)"))+
  geom_line(aes(x=time-50, y=I_t/scale1, color="Infection incidence\n(unobserved)"), linetype="dashed")+
  facet_wrap(.~tag)+
  xlab("Time")+
  scale_color_manual(values=c('black','red2'))+
  scale_y_continuous("Proportion of tests positive for influenza (ILI+)", sec.axis = sec_axis(~.*scale1, "Infection incidence"))+
  theme_bw()
