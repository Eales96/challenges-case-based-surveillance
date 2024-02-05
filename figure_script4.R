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
for(i in seq_len(8)){
  S_t <- rep((i-1)*0.05, nrow(test))
  print(i)
  
  ts_tmp<-calculate_time_series(I_t, S_t, Ts_t, Ta_t,
                             f_tau, p_tau, epsilon)
  ts_tmp$tag <- (i-1)*0.05
  ts <- rbind(ts, ts_tmp)
  
}

ts <-ts[ts$time>=50 & ts$time<250,]


###########################################################################################
#Sub-figure 1
###########################################################################################




F_tau<- f_tau
P_tau <- P_tau_calc(seq(0,30,1), beta1=1.51, beta2=2.19, beta3 = -1.1, C=3.18)

conv1 <- data.frame(tau = seq(0,30, 1),
                    F_tau1 = f_tau,
                    P_tau = P_tau)

ggplot(conv1)+
  geom_line(aes(x=tau, y=F_tau, color='F2(tau)'))+
  geom_line(aes(x=tau, y=P_tau, color='P(tau)'))+
  xlab(expression("Time since infection ("~tau~")"))+
  ylab("Probability")+
  theme_bw(base_size = 8)+
  coord_cartesian(xlim=c(0,30), ylim=c(0,1))+
  scale_color_manual("Distribution",values = c('red2', 'blue2'), labels = c(expression(F(tau)),expression(P(tau))))+
  theme(legend.position=c(0.77,0.6),
        legend.backgroun=element_rect(color='black'),
        plot.background = element_rect(color='black'))


#ggsave('figures/figure4.pdf', width=12, height=5, units= 'in')
###########################################################################################
#Sub-figure 2
###########################################################################################


plt1m<-ggplot(ts)+
  geom_line(aes(x=time, y=cs_t, color=factor(tag)))+
  geom_line(aes(x=time, y=I_t*5), linetype='dashed', color='red3')+
  scale_color_viridis_d("S(t)=", option='viridis')+
  scale_y_continuous(expression("Symptomatic cases "~C[s](t)~" (factor of "~T[s]~")"), sec.axis = sec_axis(~.*0.2, "Incidence, I(t)"))+
  xlab("Time")+
  theme_bw()+
  geom_label(aes(x = 150, y = 0.13, label = "Symptomatic cases"), fill = "white")+
  theme(legend.position="none")


plt2m<-ggplot(ts)+
  geom_line(aes(x=time, y=ca_t, color=factor(tag)))+
  geom_line(aes(x=time, y=I_t*5), linetype='dashed', color='red3')+
  scale_color_viridis_d("S(t)=", option='viridis')+
  scale_y_continuous(expression("Asymptomatic cases "~C[a](t)~" (factor of "~T[a]~")"), sec.axis = sec_axis(~.*0.2, "Incidence, I(t)"))+
  xlab("Time")+
  theme_bw()+
  geom_label(aes(x = 150, y = 0.13, label = "Asymptomatic cases"), fill = "white")+
  theme(legend.position=c(0.85,0.75),
        legend.background = element_rect(color='black'))



plt3m<-ggplot(ts)+
  geom_line(aes(x=time, y=ns_t, color=factor(tag)))+
  geom_line(aes(x=time, y=I_t*5), linetype='dashed', color='red3')+
  scale_color_viridis_d("S(t)=", option='viridis')+
  scale_y_continuous(expression("Proportion symptomatic "~n[s](t) ), sec.axis = sec_axis(~.*0.2, "Incidence, I(t)"))+
  xlab("Time")+
  theme_bw()+
  coord_cartesian(ylim=c(0,0.4))+
  geom_label(aes(x = 150, y = 0.38, label = "Proportion of population\nexhiting symptoms"), fill = "white")+
  theme(legend.position="none")
###########################################################################################



ts_a<-ts
ts_b<-ts
ts_c<-ts

ts_a$cm_t <- (ts_a$cs_t + ts_a$ca_t)
ts_a$rat <- "a"
ts_a<-ts_a[ts_a$S_t==0,]


ts_b$cm_t <- (ts_b$cs_t + 5*ts_b$ca_t)/3
ts_b$rat <- "b"

ts_c$cm_t <- (ts_c$cs_t + 0.2*ts_c$ca_t)/0.6
ts_c$rat <- "c"

ts_abc <- rbind(ts_a, ts_b, ts_c)

lin_type_vals <- c(expression(T[a](t)~"="~T[s](t)), expression(T[a](t)~"="~5*T[s](t)), expression(T[a](t)~"="~0.2*T[s](t)))



plt4m<-ggplot(ts_abc)+
  geom_line(aes(x=time, y=cm_t, color=factor(tag), linetype=rat))+
  geom_line(aes(x=time, y=I_t*5), linetype='dashed', color='red3')+
  scale_color_viridis_d("S(t)=", option='viridis')+
  scale_y_continuous(expression("All cases "~C(t)~" (factor of average testing rate)"), sec.axis = sec_axis(~./5, "Incidence, I(t)"))+
  xlab("Time")+
  theme_bw()+
  coord_cartesian(ylim=c(0.00,0.21))+
  scale_linetype_manual("Relative testing rates", labels=lin_type_vals, values=c(1,3,2))+
  geom_label(aes(x = 150, y = 0.205, label = "Cases identified through\nmass testing"), fill = "white")+
  theme(legend.position=c(0.66,0.07),
        legend.background=element_rect(color='black'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))+
  guides(color="none",
         linetype=guide_legend(nrow=1))




grid1 <- plot_grid(plt1m, plt2m, nrow=1)
grid2 <- plot_grid(plt3m, plt4m, nrow=1)

final_grid <- plot_grid(grid1, grid2, nrow=2, rel_heights = c(1,1))

ggsave('figures/figure4_alt.pdf', width=12, height=12, units= 'in')



################################################################################################
# Supplementary
###########################################################################################


ts_dif1 <- ts
ts_dif1$mult <- ts_dif1$cs_t/ts[ts$tag==0,]$cs_t
ts_dif1$add <- ts_dif1$cs_t-ts[ts$tag==0,]$cs_t

ts_dif2 <- ts
ts_dif2$mult <- ts_dif2$ca_t/ts[ts$tag==0,]$ca_t
ts_dif2$add <- ts_dif2$ca_t-ts[ts$tag==0,]$ca_t

ts_dif3 <- ts
ts_dif3$mult <- ts_dif3$ns_t/ts[ts$tag==0,]$ns_t
ts_dif3$add <- ts_dif3$ns_t-ts[ts$tag==0,]$ns_t

#ts_dif1 <- ts
#ts_dif1$mult <- ts_dif1$ca_t/ts[ts$tag==0,]$ca_t
#ts_dif1$add <- ts_dif1$ca_t-ts[ts$tag==0,]$ca_t

ts_dif1$name <- "Symptomatic cases"
ts_dif2$name <- "Asymptomatic cases"
ts_dif3$name <- "Proportion symptomatic"

ts_dif <- rbind(ts_dif1, ts_dif2, ts_dif3)


ggplot(ts_dif)+
  geom_line(aes(x=time, y=add, color=factor(tag)))+
  scale_color_viridis_d("S(t)=",option = "viridis")+
  scale_y_continuous("Additive difference")+
  xlab("Time")+
  facet_wrap(.~name,nrow=2, scales = "free_y")+
  theme_bw()+
  theme(legend.position=c(0.8,0.8),
        plot.background = element_rect(color='black'),
        legend.background = element_rect(color='black'),
        panel.grid = element_blank())


ggplot(ts_dif)+
  geom_line(aes(x=time, y=mult, color=factor(tag)))+
  scale_color_viridis_d("S(t)=",option = "viridis")+
  scale_y_continuous("Multiplicative difference")+
  xlab("Time")+
  facet_wrap(.~name,nrow=2, scales = "free_y")+
  theme_bw()+
  theme(legend.position=c(0.8,0.8),
        plot.background = element_rect(color='black'),
        legend.background = element_rect(color='black'),
        panel.grid = element_blank())


################################################################################################
# Sub-figure 4
###########################################################################################
ts_a<-ts
ts_b<-ts
ts_c<-ts

ts_a$cm_t <- (ts_a$cs_t + ts_a$ca_t)
ts_a$rat <- "a"


ts_b$cm_t <- (ts_b$cs_t + 5*ts_b$ca_t)/3
ts_b$rat <- "b"

ts_c$cm_t <- (ts_c$cs_t + 0.2*ts_c$ca_t)/0.6
ts_c$rat <- "c"

ts_abc <- rbind(ts_a, ts_b, ts_c)
ts_dif4 <- ts_abc

div<-c(ts_abc[ts_abc$tag==0&ts_abc$rat=="a",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="a",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="a",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="a",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="a",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="a",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="a",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="a",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="b",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="b",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="b",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="b",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="b",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="b",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="b",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="b",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="c",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="c",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="c",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="c",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="c",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="c",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="c",]$cm_t,
       ts_abc[ts_abc$tag==0&ts_abc$rat=="c",]$cm_t)


ts_dif4$mult <- ts_dif4$cm_t/div
ts_dif4$add <- ts_dif4$cm_t-div

ts_dif$rat='a'
ts_dif$cm_t = 0

ts_dif4$name <- "All cases"
ts_dif <- rbind(ts_dif,ts_dif4)



ts_dif$name <- factor(ts_dif$name, levels=c("Symptomatic cases", "Asymptomatic cases", "Proportion symptomatic", "All cases"))

additive<-ggplot(ts_dif)+
  geom_line(data=ts_dif[ts_dif$rat=="a",], aes(x=time, y=add, color=factor(tag)))+
  geom_line(data=ts_dif[ts_dif$rat=="a" &ts_dif$S_t==0,], aes(x=time, y=add, color=factor(tag)))+
  geom_line(data=ts_dif[ts_dif$rat %in% c("b", "c"),], aes(x=time, y=add, color=factor(tag), linetype=rat))+
  scale_color_viridis_d("S(t)=",option = "viridis")+
  scale_y_continuous("Additive difference")+
  scale_linetype_manual("Relative testing rates",labels=lin_type_vals[2:3], values=c(3,2))+
  xlab("Time")+
  facet_wrap(.~name,nrow=2, scales = "free_y")+
  theme_bw()+
  theme(legend.position="bottom",
        plot.background = element_rect(color='black'),
        legend.background = element_rect(color='black'),
        panel.grid = element_blank())+
  guides(color=guide_legend(nrow=2))





ggsave('figures/figure4_altsup1.pdf', width=8, height=8, units= 'in')



multi<-ggplot(ts_dif)+
  geom_line(data=ts_dif[ts_dif$rat=="a",], aes(x=time, y=mult, color=factor(tag)))+
  geom_line(data=ts_dif[ts_dif$rat=="a" &ts_dif$S_t==0,], aes(x=time, y=mult, color=factor(tag)))+
  geom_line(data=ts_dif[ts_dif$rat %in% c("b", "c"),], aes(x=time, y=mult, color=factor(tag), linetype=rat))+
  scale_color_viridis_d("S(t)=",option = "viridis")+
  scale_y_continuous("Multiplicative difference")+
  scale_linetype_manual("Relative testing rates",labels=lin_type_vals[2:3], values=c(3,2))+
  xlab("Time")+
  facet_wrap(.~name,nrow=2, scales = "free_y")+
  theme_bw()+
  theme(legend.position="bottom",
        plot.background = element_rect(color='black'),
        legend.background = element_rect(color='black'),
        panel.grid = element_blank())+
  guides(color=guide_legend(nrow=2))





ggsave('figures/figure4_altsup2.pdf', width=8, height=8, units= 'in')



#################################################################################################




###################################################################################################################################

