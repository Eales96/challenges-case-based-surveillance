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


s <- seq(0,1,0.1)

as <- 1 - s

df_as <- data.frame(s = s,
                    rel = 1-s,
                    tag = 0)




plt1m <- ggplot(df_as, aes(x=s, y=rel))+
  geom_line(color='blue')+
  theme_bw()+
  xlab("Background proportion symptomatic S(t)")+
  ylab("Asymptomatic cases (relative to value at S(t)=0)")+
  coord_cartesian(xlim=c(0.009,0.195), ylim=c(0.7,1))+
  geom_hline(yintercept = 1)+
  theme(legend.position = "none",
        legend.background = element_rect(color='black'))+
  geom_label(aes(x=0.05, y=0.8, label="Gradient=-1"), color='black')



plt1 <- ggdraw()+
  draw_plot(plt1m)+
  labs(tag = "A")+
  theme(plot.tag.position = c(0.02,0.98))



###########################################################################################
#Sub-figure 2
###########################################################################################

s <- seq(0, 0.3, 0.01)
cs <- 10**seq(-1, 1, 0.1)

df_s<- expand.grid(s, cs)

colnames(df_s) <- c("s", "tag")
df_s$rel <- 1+df_s$s*df_s$tag

leg_tit <- expression(frac(integral(I(t-tau)*(1-F(tau))*P(tau)*d~tau, 0, infinity), integral(I(t-tau)*F(tau)*P(tau)*d~tau, 0, infinity)))
grad_tit <- expression("Gradient="~frac(integral(I(t-tau)*(1-F(tau))*P(tau)*d~tau, 0, infinity), integral(I(t-tau)*F(tau)*P(tau)*d~tau, 0, infinity)))


plt2m<-ggplot(df_s, aes(x=s, y=rel, col=tag, group=factor(tag)))+
  geom_line()+
  theme_bw()+
  scale_color_viridis_c(leg_tit, option="turbo", trans="log10")+
  geom_vline(xintercept = 0.02, linetype="dashed")+
  xlab("Background proportion symptomatic S(t)")+
  ylab("Symptomatic cases (relative to value at S(t)=0)")+
  coord_cartesian(xlim=c(0.009,0.195))+
  geom_hline(yintercept = 1)+
  theme(legend.position = "none",
        legend.background = element_rect(color='black'))+
  guides(color = guide_colourbar(label.position = "right"))+
  geom_label(aes(x=0.153, y=3.75, label=paste(grad_tit)), parse=T, color='black')



plt2i<-ggplot(df_s[df_s$s==0.02,], aes(x=tag, y=rel))+
  geom_line(linetype="dashed")+
  geom_point(aes(x=tag, y=rel, col=tag))+
  ylab("Relative value")+
  xlab(leg_tit)+
  theme_bw(base_size = 10)+
  scale_color_viridis_c(option="turbo", trans="log10")+
  scale_x_log10()+
  theme(legend.position = "none",
        plot.background = element_rect(color="black"))



plt2 <- ggdraw()+
  draw_plot(plt2m)+
  draw_plot(plt2i, x=0.18, y=0.58, width=0.36, height=0.4)+
  ggtitle("Test")+
  labs(tag = "B")+
  theme(plot.tag.position = c(0.02,0.98))


###########################################################################################
#Sub-figure 3
###########################################################################################

s <- seq(0, 0.3, 0.01)
ns <- 10**seq(-3, -1, 0.1)

df_ns<- expand.grid(s, ns)



colnames(df_ns) <- c("s", "tag")
df_ns$rel <- 1 + df_ns$s + (df_ns$s/df_ns$tag)


leg_tit <- expression(integral(I(t-tau)*F(tau)*d~tau, 0, infinity))
grad_tit <- expression("Gradient="~frac(1-integral(I(t-tau)*F(tau)*d~tau, 0, infinity),integral(I(t-tau)*F(tau)*d~tau, 0, infinity)))

plt3m<-ggplot(df_ns, aes(x=s, y=rel, col=tag, group=factor(tag)))+
  geom_line()+
  theme_bw()+
  scale_color_viridis_c(leg_tit, option="turbo", trans="log10")+
  geom_vline(xintercept = 0.02, linetype="dashed")+
  xlab("Background proportion symptomatic S(t)")+
  ylab("Proportion symptomatic (relative to value at S(t)=0)")+
  coord_cartesian(xlim=c(0.009,0.195))+
  geom_hline(yintercept = 1)+
  theme(legend.position = "none",
        legend.background = element_rect(color='black'))+
  guides(color = guide_colourbar(label.position = "right"))+
  geom_label(aes(x=0.15, y=280, label=paste(grad_tit)), parse=T, color='black')


plt3i<-ggplot(df_ns[df_ns$s==0.02,], aes(x=tag, y=rel))+
  geom_line(linetype="dashed")+
  geom_point(aes(x=tag, y=rel, col=tag))+
  ylab("Relative value")+
  xlab(leg_tit)+
  theme_bw(base_size = 10)+
  scale_color_viridis_c(option="turbo", trans="log10")+
  scale_x_log10()+
  theme(legend.position = "none",
        plot.background = element_rect(color="black"))



plt3 <- ggdraw()+
  draw_plot(plt3m)+
  draw_plot(plt3i, x=0.2, y=0.58, width=0.36, height=0.4)+
  labs(tag = "C")+
  theme(plot.tag.position = c(0.02,0.98))


###########################################################################################
#Sub-figure 3
###########################################################################################

s <- seq(0, 0.3, 0.01)
mt_a <- 10**seq(-1, 1, 0.2)
mt_b <- c(5,0.2)

df_mt<- expand.grid(s, mt_a, mt_b)


colnames(df_mt) <- c("s", "tag_a", "tag_b")
df_mt$rel <- 1 + (df_mt$s*df_mt$tag_a*(1-df_mt$tag_b)/(1+df_mt$tag_a*df_mt$tag_b))



grad_tit <- expression("Gradient="~
                         frac( (1-frac(T[a],T[s])) *frac(integral(I(t-tau)*(1-F(tau))*P(tau)*d~tau, 0, infinity), integral(I(t-tau)*F(tau)*P(tau)*d~tau, 0, infinity)),
                               1 + (frac(T[a],T[s])) *frac( integral(I(t-tau)*(1-F(tau))*P(tau)*d~tau, 0, infinity), integral(I(t-tau)*F(tau)*P(tau)*d~tau, 0, infinity)) ))


lin_type_vals <- c(expression(T[a](t)~"="~0.2*T[s](t)), expression(T[a](t)~"="~5*T[s](t)))

plt4m<-ggplot(df_mt)+
  geom_line(aes(x=s, y=rel, col=tag_a, linetype=factor(tag_b) ,group=interaction(tag_a, tag_b)))+
  theme_bw()+
  scale_color_viridis_c(leg_tit, option="turbo", trans="log10")+
  scale_linetype_manual("Relative testing rates",labels=lin_type_vals, values=c(2,3))+
  geom_vline(xintercept = 0.02, linetype="dashed")+
  xlab("Background proportion symptomatic S(t)")+
  ylab("All cases (relative to value at S(t)=0)")+
  coord_cartesian(xlim=c(0.009,0.195), ylim=c(0.1,1.6))+
  geom_hline(yintercept = 1)+
  theme(legend.position = c(0.845,0.2),
        legend.background = element_rect(color='black'))+
  geom_label(aes(x=0.07, y=1.5, label=paste(grad_tit)), parse=T)+
  guides(color = "none")


leg_tit <- expression(frac(integral(I(t-tau)*(1-F(tau))*P(tau)*d~tau, 0, infinity), integral(I(t-tau)*F(tau)*P(tau)*d~tau, 0, infinity)))

plt4i<-ggplot(df_mt[df_mt$s==0.02,])+
  geom_line(aes(x=tag_a, y=rel,group=factor(tag_b), linetype=factor(tag_b)))+
  geom_point(aes(x=tag_a, y=rel, col=tag_a, group=factor(tag_b)))+
  ylab("Relative value")+
  xlab(leg_tit)+
  theme_bw(base_size = 10)+
  scale_linetype_manual("Relative testing rates",labels=lin_type_vals, values=c(2,3))+
  scale_color_viridis_c(option="turbo", trans="log10")+
  scale_x_log10()+
  theme(legend.position = "none",
        legend.background = element_rect(color='black'),
        plot.background = element_rect(color="black"))+
  guides(color = "none")


plt4 <- ggdraw()+
  draw_plot(plt4m)+
  draw_plot(plt4i, x=0.2, y=0.1, width=0.4, height=0.45)+
  labs(tag = "D")+
  theme(plot.tag.position = c(0.02,0.98))



#################################################################################
grd_plt <- plot_grid(plt1, plt2, plt3, plt4, nrow = 2)

ggsave('figures/figure5.pdf', width=12, height=12, units= 'in')

ts1<-c(rep(0,1600),rep(1,26))
ts2<-c(rep(0,1600),rep(1,13))
wilcox.test(ts1, ts2)
