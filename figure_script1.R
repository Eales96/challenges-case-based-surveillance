setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/BiasesInCaseData")

source('R/functions.R')
library(ggplot2)



# Sub figure 1 (Symptom onset distributions)
symp_onset <- data.frame()

tau <- seq(0,10,0.1)
tag1 <- rep("2 days",101)
tag2 <- rep("4 days",101)

symp_onset <- data.frame(tau=c(tau,tau))
symp_onset$tag <- c(tag1, tag2)
symp_onset$pdf <- c(0.8*dgamma(tau, shape=8, rate=4), 0.8*dgamma(tau, shape=4, rate=1))
symp_onset$cdf <- c(0.8*pgamma(tau, shape=8, rate=4), 0.8*pgamma(tau, shape=4, rate=1))


plt1m <- ggplot(data=symp_onset)+
  geom_line(aes(x=tau, y= pdf, color=tag))+
  ylab(expression("Probability density function, "~f[1](tau)))+
  xlab("Time to symptom onset")+
  #geom_hline(yintercept = 0.4, linetype='dashed')+
  scale_color_brewer("Mean time to\nsymptom onset",palette = "Dark2")+
  theme_bw()+
  theme(legend.position = c(0.8,0.3),
        legend.background = element_rect(color="black"),
        legend.title =element_text(size=8),
        legend.text =element_text(size=8))


plt1i <- ggplot(data=symp_onset)+
  geom_line(aes(x=tau, y=cdf, color=tag))+
  ylab("Cumulative density\nfunction")+
  xlab("Time to symptom onset")+
  geom_hline(yintercept = 0.8, linetype='dashed', color='blue3')+
  scale_color_brewer("Mean time\nto symptom onset",palette = "Dark2")+
  theme_bw(base_size = 8)+
  coord_cartesian(ylim = c(0,1))+
  theme(legend.position = "none",
        plot.background = element_rect(color = 'black'))


library(cowplot)

plt1 <- ggdraw()+
  draw_plot(plt1m)+
  draw_plot(plt1i, x=0.5, y=0.55, width=0.45, height=0.4)+
  labs(tag = "A")+
  theme(plot.tag.position = c(0.02,0.98))




##############################################################################################################
# Sub figure 2 (Symptom duration distribution)
symp_dur <- data.frame()

tau <- seq(0,10,0.1)
tag1 <- rep("2 days",101)
tag2 <- rep("5 days",101)

symp_dur <- data.frame(tau=c(tau,tau))
symp_dur$tag <- c(tag1, tag2)
symp_dur$pdf <- c(dgamma(tau, shape=20, rate=10), dgamma(tau, shape=50, rate=10))
symp_dur$cdf <- c(pgamma(tau, shape=20, rate=10), pgamma(tau, shape=50, rate=10))



plt2m<-ggplot(data=symp_dur)+
  geom_line(aes(x=tau, y= 1-cdf, linetype = tag))+
  ylab(expression("Proportion still symptomatic, "~f[2](theta)))+
  xlab("Time since symptom onset")+
  scale_linetype_manual("Mean symptom duration", values=c(3,2))+
  theme_bw()+
  theme(legend.position = c(0.78,0.32),
        legend.background = element_rect(color="black"),
        legend.title =element_text(size=8),
        legend.text =element_text(size=8)) 

plt2i<-ggplot(data=symp_dur)+
  geom_line(aes(x=tau, y= pdf, linetype = tag))+
  ylab("Probability density\nfunction")+
  xlab("Symptom duration")+
  scale_linetype_manual("Symptom duration", values=c(3,2))+
  theme_bw(base_size = 8)+
  theme(legend.position = "none",
        plot.background = element_rect(color='black'))



plt2 <- ggdraw()+
  draw_plot(plt2m)+
  draw_plot(plt2i, x=0.57, y=0.57, width=0.4, height=0.4)+
  labs(tag = "B")+
  theme(plot.tag.position = c(0.02,0.98))





################################################################################################################
wt_sum <- function(pdf1, pdf2){
  dist <- rep(0, length(pdf1)+length(pdf2)-1)
  for(i in seq_len(length(pdf1))){
    dist[i:(i-1+length(pdf2))] <- dist[i:(i-1+length(pdf2))] + pdf1[i]*pdf2
  }
  dist*0.8/sum(pdf1)
}

##################################################################################################################
#Sub-figure 3 (Convolutions)

df1 <- data.frame(conv=wt_sum(symp_onset[symp_onset$tag=="2 days",]$pdf,
                              1-symp_dur[symp_dur$tag=="2 days",]$cdf),
                  tau = seq(0,20,0.1),
                  onset = "2 days",
                  dur = "2 days")


df2 <- data.frame(conv=wt_sum(symp_onset[symp_onset$tag=="2 days",]$pdf,
                              1-symp_dur[symp_dur$tag=="5 days",]$cdf),
                  tau = seq(0,20,0.1),
                  onset = "2 days",
                  dur = "5 days")


df3 <- data.frame(conv=wt_sum(symp_onset[symp_onset$tag=="4 days",]$pdf,
                              1-symp_dur[symp_dur$tag=="2 days",]$cdf),
                  tau = seq(0,20,0.1),
                  onset = "4 days",
                  dur = "2 days")


df4 <- data.frame(conv=wt_sum(symp_onset[symp_onset$tag=="4 days",]$pdf,
                              1-symp_dur[symp_dur$tag=="5 days",]$cdf),
                  tau = seq(0,20,0.1),
                  onset = "4 days",
                  dur = "5 days")



df <- rbind(df1, df2, df3, df4)
plot(df$conv)
plt3<-ggplot(data=df)+
  geom_line(aes(x=tau, y= conv, linetype=dur, color=onset))+
  xlab("Time since infection")+
  ylab(expression("Proportion exhibiting symptoms, "~F(tau)))+
  scale_linetype_manual("Mean symptom duration", values=c(3,2))+
  scale_color_brewer("Mean time to\nsymptom onset",palette = "Dark2")+
  theme_bw()+
  labs(tag="C")+
  geom_hline(yintercept = 0.8, linetype='dashed', color='blue3')+
  theme(legend.position = c(0.76,0.5),
        legend.background=element_rect(color='black'),
        plot.tag.position = c(0.02,0.98),
        legend.title =element_text(size=8),
        legend.text =element_text(size=8))



plot_grid(plt1, plt2, plt3, nrow=1, rel_widths = c(1,1,1))
ggsave('figures/figure1.pdf', width=12, height=4, units= 'in')
