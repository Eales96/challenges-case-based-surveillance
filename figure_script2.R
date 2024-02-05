#setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/BiasesInCaseData")

source('R/functions.R')
library(ggplot2)



# Sub figure 1 (Symptom onset distributions)
symp_onset <- data.frame()

tau <- seq(0,30,1)

P_tau <- function(tau, beta1, beta2, beta3, C){
  x = tau - C
  val = beta1 + x*beta2 +x*beta2*beta3*(x>0)
  exp(val)/(1+exp(val))
}

df1<-data.frame(tau= tau,
                pos =P_tau(tau, beta1=1.51, beta2=2.19, beta3 = -1.1, C=2.18),
                tag = "a(tau)")

df2<-data.frame(tau= tau,
                pos =0.05,
                tag = "Epsilon")

df3<- data.frame(tau=seq(32,38,1),
                 pos=0.05,
                 tag="Epsilon")
df <- rbind(df1, df2)


mn1<-mean(df1[df1$tau>0 &df1$tau<6,]$pos+df2[df2$tau>0 &df2$tau<6,]$pos)
mn2<-mean(df1[df1$tau>5 &df1$tau<11,]$pos+df2[df2$tau>5 &df2$tau<11,]$pos)
mn3<-mean(df1[df1$tau>10 &df1$tau<26,]$pos+df2[df2$tau>10 &df2$tau<26,]$pos)
mn<-data.frame(sensitivity=c(mn1, mn2, mn3),
           period=c("1","2","3"))

plt1i<-ggplot(data=mn, aes(x=period, y=sensitivity,color=period))+
  geom_point(size=4)+
  scale_color_manual(values = c('red2','blue2', 'orange2'))+
  theme_bw(base_size = 8)+
  xlab("Example sampling period")+
  ylab("Effective sensitivity")+
  coord_cartesian(ylim=c(0,1))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.background = element_rect(color='black'))
  


plt1m<-ggplot(data=df)+
  annotate("rect", xmin=1, xmax=5, ymin=0, ymax=Inf, alpha=0.2, fill="red2")+
  annotate("text", x = 3, y = 0.9, label = "Example\nsampling\nperiod 1", size=3)+
  annotate("rect", xmin=6, xmax=10, ymin=0, ymax=Inf, alpha=0.2, fill="blue2")+
  annotate("text", x = 8, y = 0.9, label = "Example\nsampling\nperiod 2", size=3)+
  annotate("rect", xmin=11, xmax=25, ymin=0, ymax=Inf, alpha=0.2, fill="orange2")+
  annotate("text", x = 18, y = 0.9, label = "Example\nsampling\nperiod 3",size=3)+
  geom_area(aes(x=tau, y=pos, fill=tag), alpha=0.8)+
  geom_area(data=df3,aes(x=tau, y=pos, fill=tag))+
  geom_line(data=df1,aes(x=tau, y=pos+0.05), color='black')+
  geom_line(data=df3,aes(x=tau, y=pos), color='black')+
  ylab("Probability of testing positive")+
  xlab("Time since infection (days)")+
  scale_fill_brewer("Component",palette = "Dark2", labels=parse(text=c(paste("   P(tau)"), paste("epsilon"))))+
  theme_bw()+
  scale_x_continuous(breaks=c(0,10,20,30, 30.5, 31,31.5,32.,32.5,33,33.5,34,34.5,35),labels=c(0,10,20,30,"",".....","",parse(text=paste("infinity")),"","","","","",""))+
  coord_cartesian(xlim=c(1.5,33))+
  theme(legend.position = c(0.9,0.5),
        legend.background = element_rect(color="black"))






library(cowplot)

plt1 <- ggdraw()+
  draw_plot(plt1m)+
  draw_plot(plt1i, x=0.62, y=0.62, width=0.35, height=0.35)
plt1

ggsave('figures/figure2.pdf', width=7, height=6, units= 'in')


