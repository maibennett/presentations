rm(list = ls())
cat("\014")

##########################################################
######## GIF RDD  for Template Matching ################## 
##########################################################

#devtools::install_github("dgrtwo/gganimate")

library(tidyverse)
library(gganimate)
library(ggthemes)
library(gridExtra)
library(cowplot)

set.seed(100)

df_pre <- data.frame(xaxisTime=runif(300)*20) %>%
  mutate(Y = .3*xaxisTime+rnorm(300),
         state="1",
         groupX=floor(xaxisTime)+.5,
         groupLine=floor(xaxisTime),
         cutLine=rep(c(9,11),150)) %>%
  group_by(groupX) %>%
  mutate(mean_Y=mean(Y)) %>%
  ungroup() %>%
  arrange(groupX)

df_pre$id = seq(1,nrow(df_pre),1)

df_post <- data.frame(xaxisTime=runif(300)*20) %>%
  mutate(Y = .3*xaxisTime+5*(xaxisTime>10)-.15*xaxisTime*(xaxisTime>10)+rnorm(300),
         state="1",
         groupX=floor(xaxisTime)+.5,
         groupLine=floor(xaxisTime),
         cutLine=rep(c(9,11),150)) %>%
  group_by(groupX) %>%
  mutate(mean_Y=mean(Y)) %>%
  ungroup() %>%
  arrange(groupX)

y1 = df_pre$Y + exp(-2*(df_pre$xaxisTime-3)) + 1 - 0.04*df_pre$xaxisTime
poly_fit_y1 <- as.data.frame(lowess(df_pre$xaxisTime, y1, f=0.4))

df_post2 <- data.frame(xaxisTime=runif(300)*20) %>%
  mutate(Y = .3*xaxisTime+(xaxisTime<10)*(poly_fit_y1$y-df_pre$Y)+rnorm(300),
         state="1",
         groupX=floor(xaxisTime)+.5,
         groupLine=floor(xaxisTime),
         cutLine=rep(c(9,11),150)) %>%
  group_by(groupX) %>%
  mutate(mean_Y=mean(Y)) %>%
  ungroup() %>%
  arrange(groupX)

df_post2_detrend <- data.frame(xaxisTime=runif(300)*20) %>%
  mutate(Y = rnorm(300)+2.8 + (xaxisTime<10)*(poly_fit_y1$y-df_pre$Y),
         state="1",
         groupX=floor(xaxisTime)+.5,
         groupLine=floor(xaxisTime),
         cutLine=rep(c(9,11),150)) %>%
  group_by(groupX) %>%
  mutate(mean_Y=mean(Y)) %>%
  ungroup() %>%
  arrange(groupX)

df_post$id = seq(1,nrow(df_post),1)
df_post2$id = seq(1,nrow(df_post2),1)
df_post2_detrend$id = seq(1,nrow(df_post2_detrend),1)

dir_data = dirname(rstudioapi::getSourceEditorContext()$path)
save(df_pre,df_post, df_post2, df_post2_detrend, file=paste0(dir_data,"/data_plot.Rdata"))


##########################################################
rm(list = ls())
cat("\014")

##########################################################
######## GIF RDD  for Template Matching ################## 
##########################################################

#devtools::install_github("dgrtwo/gganimate")

library(tidyverse)
library(gganimate)
library(ggthemes)
library(gridExtra)
library(pBrackets)

dir_data = dirname(rstudioapi::getSourceEditorContext()$path)
load(paste0(dir_data,"/data_plot.Rdata"))

par(mfcol=c(1,1),xpd=FALSE)

poly_fit_y0 <- as.data.frame(lowess(df_pre$xaxisTime, df_pre$Y, f=0.4))

y1 = df_pre$Y + exp(-2*(df_pre$xaxisTime-3)) + 1 - 0.04*df_pre$xaxisTime
poly_fit_y1 <- as.data.frame(lowess(df_pre$xaxisTime, y1, f=0.4))

#Pre-period
g1 <- ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_line(data = poly_fit_y0, aes(x=x, y=y, color='deepskyblue3'), size=1)+
  geom_line(data = poly_fit_y1, aes(x=x, y=y, color = 'darkorchid3'), size=1, lty=2)+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_colour_identity(name='Potential Outcomes',breaks=c('deepskyblue3','darkorchid3'),
                        labels = c(bquote(Y^{(0)}~(R)),bquote(Y^{(1)}~(R))), guide="legend") +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),legend.position="bottom")+
  ylim(-2,10)


#Post-period
g2 <- ggplot(df_post2,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_line(data = poly_fit_y0[poly_fit_y0$x>10,], aes(x=x, y=y,color='deepskyblue3'), size=1)+
  geom_line(color='deepskyblue3',data = poly_fit_y0[poly_fit_y0$x<10,], aes(x=x, y=y), size=1,lty=2)+
  geom_line(data = poly_fit_y1[poly_fit_y1$x>10,], aes(x=x, y=y, color='darkorchid3'), size=1, lty=2)+
  geom_line(color='darkorchid3',data = poly_fit_y1[poly_fit_y1$x<10,], aes(x=x, y=y), size=1, lty=1)+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Post-intervention")+
  scale_colour_identity(name='Potential Outcomes',breaks=c('deepskyblue3','darkorchid3'),
                        labels = c(bquote(Y^{(0)}~(R)),bquote(Y^{(1)}~(R))), guide="legend") +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),legend.position="bottom") +
  ylim(-2,10)

library(cowplot)
plot_grid(g1, g2)


y1 = df_post2_detrend$Y + exp(-2*(df_post2_detrend$xaxisTime-3)) + 1 - 0.04*df_post2_detrend$xaxisTime
poly_fit_y1 <- as.data.frame(lowess(df_post2_detrend$xaxisTime, y1, f=0.4))

#Randomly choose some points for matching:
set.seed(1)
T_matched = sample(df_post2_detrend$id[df_post2_detrend$xaxisTime<10 & df_post2_detrend$xaxisTime>5],size=30) 
C_matched = sample(df_post2_detrend$id[df_post2_detrend$xaxisTime>10 & df_post2_detrend$xaxisTime<15],size=30)

#Post-period after matching
ggplot(df_post2_detrend[df_post2_detrend$xaxisTime>=5 & df_post2_detrend$xaxisTime<=15,],
       aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="light grey", color="light grey")+
  geom_point(data=df_post2_detrend[which(df_post2_detrend$id %in% T_matched),],
             aes(x=xaxisTime,y=Y,color="springgreen4"),pch=1) + 
  geom_point(data=df_post2_detrend[which(df_post2_detrend$id %in% C_matched),],
             aes(x=xaxisTime,y=Y,color="springgreen4"),pch=2) + 
  geom_line(data = poly_fit_y0[poly_fit_y0$x>10,], aes(x=x, y=2.8,color='deepskyblue3'), size=1)+
  geom_line(color='deepskyblue3',data = poly_fit_y0[poly_fit_y0$x<10,], aes(x=x, y=2.8), size=1,lty=2)+
  geom_line(data = poly_fit_y1[poly_fit_y1$x>10,], aes(x=x, y=y, color='darkorchid3'), size=1, lty=2)+
  geom_line(color='darkorchid3',data = poly_fit_y1[poly_fit_y1$x<10,], aes(x=x, y=y), size=1, lty=1)+
  geom_vline(aes(xintercept=10),linetype='dashed', color="black") +
  geom_vline(aes(xintercept=5),linetype='dashed',color='springgreen4',lwd=1) +
  geom_vline(aes(xintercept=15),linetype='dashed',color='springgreen4',lwd=1) + theme_bw()+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Post-intervention after Matching")+
  scale_colour_identity(name='Potential Outcomes',breaks=c('deepskyblue3','darkorchid3',"springgreen4"),
                        labels = c("Y0|X","Y1|X","Matched obs"),lty=c(1,1,NA),
                        pch=c(NA,NA,1), guide="legend") +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),legend.position="bottom") +
  ylim(0,10) + xlim(4,16)


lm_fit1 <- lm(Y ~ xaxisTime, data=df_post[df_post$xaxisTime<10,])
summary(lm_fit1)

lm_fit2 <- lm(Y ~ xaxisTime, data=df_post[df_post$xaxisTime>10,])
summary(lm_fit2)

predicted_df1 <- data.frame(xaxisTime = df_post$xaxisTime[df_post$xaxisTime<10],
                            pred1 = predict(lm_fit1))
predicted_df2 <- data.frame(xaxisTime = df_post$xaxisTime[df_post$xaxisTime>10],
                            pred2 = predict(lm_fit2))

ggplot(df_post,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_line(color='dark grey',data = predicted_df1, aes(x=xaxisTime, y=pred1), size=1)+
  geom_line(color='dark grey',data = predicted_df2, aes(x=xaxisTime, y=pred2), size=1)+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  xlab("Running variable") + ylab("Outcome")+
  annotate("rect", xmin = 9.2, xmax = 10.8, ymin = 10.5, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 10, y = 11, label = "Cutoff", hjust=0.5, size=5)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12))

  

p1 = c(10,lm_fit1$coefficients[1]+10*lm_fit1$coefficients[2])
p2 = c(10,lm_fit2$coefficients[1]+10*lm_fit2$coefficients[2])

ggplot(df_post,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_line(color='dark grey',data = predicted_df1, aes(x=xaxisTime, y=pred1), size=1)+
  geom_line(color='dark grey',data = predicted_df2, aes(x=xaxisTime, y=pred2), size=1)+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_segment(aes(x = p1[1], y = p1[2], xend = p2[1], yend = p2[2]), col="#eb811b", size=1.5)+
  xlab("Running variable") + ylab("Outcome")+
  annotate("rect", xmin = 9.2, xmax = 10.8, ymin = 10.5, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 10, y = 11, label = "Cutoff", hjust=0.5, size=5)+
  annotate("text", x = 10.2, y = 4, label = "LATE", hjust=0, size=5,col="#eb811b")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12))


ggplot(df_post,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_line(color='#eb811b',data = predicted_df1, aes(x=xaxisTime, y=pred1), size=1)+
  geom_line(color='#eb811b',data = predicted_df2, aes(x=xaxisTime, y=pred2), size=1)+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_segment(aes(x = p1[1], y = p1[2], xend = p2[1], yend = p2[2]), col="dark grey", size=1.5)+
  xlab("Running variable") + ylab("Outcome")+
  annotate("rect", xmin = 9.2, xmax = 10.8, ymin = 10.5, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 10, y = 11, label = "Cutoff", hjust=0.5, size=5)+
  annotate("text", x = 10.2, y = 4, label = "LATE", hjust=0, size=5,col="dark grey")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12))



ggplot(df_post,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_line(color='lightgrey',data = predicted_df1, aes(x=xaxisTime, y=pred1), size=1, linetype='dashed')+
  geom_line(color='lightgrey',data = predicted_df2, aes(x=xaxisTime, y=pred2), size=1, linetype='dashed')+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_vline(aes(xintercept=5),linetype='dashed',color='#eb811b',lwd=2)+
  geom_vline(aes(xintercept=15),linetype='dashed',color='#eb811b',lwd=2)+
  geom_segment(aes(x = p1[1], y = p1[2], xend = p2[1], yend = p2[2]), col="dark grey", size=1.5)+
  xlab("Running variable") + ylab("Outcome")+
  annotate("rect", xmin = 9.2, xmax = 10.8, ymin = 10.5, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 10, y = 11, label = "Cutoff", hjust=0.5, size=5)+
  annotate("text", x = 10.2, y = 4, label = "LATE", hjust=0, size=5,col="dark grey")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12))


###### Plot the method:

par(mfrow=c(2,1),xpd=FALSE)

lm_fit1p <- lm(Y ~ xaxisTime, data=df_pre[df_pre$xaxisTime<10,])
summary(lm_fit1p)

lm_fit2p <- lm(Y ~ xaxisTime, data=df_pre[df_pre$xaxisTime>10,])
summary(lm_fit2p)

predicted_df1p <- data.frame(xaxisTime = df_pre$xaxisTime[df_pre$xaxisTime<10],
                            pred1 = predict(lm_fit1p))
predicted_df2p <- data.frame(xaxisTime = df_pre$xaxisTime[df_pre$xaxisTime>10],
                            pred2 = predict(lm_fit2p))

plot1 = ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_line(color='dark grey',data = predicted_df1p, aes(x=xaxisTime, y=pred1), size=1)+
  geom_line(color='dark grey',data = predicted_df2p, aes(x=xaxisTime, y=pred2), size=1)+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  annotate("rect", xmin = 8.3, xmax = 11.8, ymin = 10.6, ymax = 11.4,
           color="black", fill="white") +
  annotate("text", x = 10, y = 11, label = "Cutoff", hjust=0.5, size=5)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face = "bold"))


plot2 = ggplot(df_post,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_line(color='dark grey',data = predicted_df1, aes(x=xaxisTime, y=pred1), size=1)+
  geom_line(color='dark grey',data = predicted_df2, aes(x=xaxisTime, y=pred2), size=1)+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Post-intervention")+
  annotate("rect", xmin = 8.3, xmax = 11.8, ymin = 10.6, ymax = 11.4,
           color="black", fill="white") +
  annotate("text", x = 10, y = 11, label = "Cutoff", hjust=0.5, size=5)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))

grid.arrange(plot1, plot2, nrow = 1)

### 2

plot1 = ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_line(color='dark grey',data = predicted_df1p, aes(x=xaxisTime, y=pred1), size=1)+
  geom_line(color='dark grey',data = predicted_df2p, aes(x=xaxisTime, y=pred2), size=1)+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  annotate("rect", xmin = 8.3, xmax = 11.8, ymin = 10.6, ymax = 11.4,
           color="black", fill="white") +
  annotate("text", x = 10, y = 11, label = "Cutoff", hjust=0.5, size=5)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))


plot2 = ggplot(df_post,aes(y=Y,x=xaxisTime, label=NA), fill="grey")+geom_point(fill="grey", color="grey")+
  geom_line(color='grey',data = predicted_df1, aes(x=xaxisTime, y=pred1), size=1)+
  geom_line(color='grey',data = predicted_df2, aes(x=xaxisTime, y=pred2), size=1)+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Post-intervention")+
  annotate("rect", xmin = 8.3, xmax = 11.8, ymin = 10.6, ymax = 11.4,
           color="dark grey", fill="white") +
  annotate("text", x = 10, y = 11, label = "Cutoff", hjust=0.5, size=5, col="grey")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0), colour="grey"),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0), colour="grey"),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5), 
        plot.title = element_text(colour = "grey"))

grid.arrange(plot1, plot2, nrow = 1)


### 3

ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_vline(aes(xintercept=9),linetype='dashed', colour="dark grey", size=1) +
  geom_vline(aes(xintercept=11),linetype='dashed', colour="dark grey", size=1) + 
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "BW1", hjust=0.5, size=4, col="black")+
  geom_segment(aes(x = 9.4, y = 11, xend = 9, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 11, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))



### 4

ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_vline(aes(xintercept=9),linetype='dashed', colour="dark grey", size=1) +
  geom_vline(aes(xintercept=11),linetype='dashed', colour="dark grey", size=1) + 
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "BW1", hjust=0.5, size=4, col="black")+
  geom_segment(aes(x = 9.4, y = 11, xend = 9, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 11, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="darkorchid3", shape=22, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="deepskyblue3", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="gold2", shape=24, size=3)



### 5

ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_vline(aes(xintercept=9),linetype='dashed', colour="dark grey", size=1) +
  geom_vline(aes(xintercept=11),linetype='dashed', colour="dark grey", size=1) + 
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "BW1", hjust=0.5, size=4, col="black")+
  geom_segment(aes(x = 9.4, y = 11, xend = 9, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 11, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="darkorchid3", shape=22, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="deepskyblue3", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="gold2", shape=24, size=3)


### 6

ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_vline(aes(xintercept=9),linetype='dashed', colour="dark grey", size=1) +
  geom_vline(aes(xintercept=11),linetype='dashed', colour="dark grey", size=1) + 
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "BW1", hjust=0.5, size=4, col="black")+
  geom_segment(aes(x = 9.4, y = 11, xend = 9, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 11, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="darkorchid3", shape=22, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="deepskyblue3", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="gold2", shape=24, size=3)+
  #matched points
  geom_point(aes(y=c(df_pre$Y[df_pre$xaxisTime<10],rep(NA,sum(df_pre$xaxisTime>=10))),
                 x=c(df_pre$xaxisTime[df_pre$xaxisTime<10],rep(NA,sum(df_pre$xaxisTime>=10)))), 
             color="light grey", fill="light grey")+
  geom_point(aes(y=3.1,x=9.1),
             color="black", fill="darkorchid3", shape=22, size=2.5)+
  geom_point(aes(y=3.5,x=9.4),
             color="black", fill="deepskyblue3", shape=23, size=2.5)+
  geom_point(aes(y=2.6,x=9.7),
             color="black", fill="gold2", shape=24, size=2.5)




### 7

ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="light grey")+geom_point(fill="light grey", color="light grey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_vline(aes(xintercept=9),linetype='dashed', colour="dark grey", size=1) +
  geom_vline(aes(xintercept=11),linetype='dashed', colour="dark grey", size=1) + 
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "BW1", hjust=0.5, size=4, col="black")+
  geom_segment(aes(x = 9.4, y = 11, xend = 9, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 11, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="darkorchid3", shape=22, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="deepskyblue3", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="gold2", shape=24, size=3)+
  #matched points control
  geom_point(aes(y=3.1,x=9.1),
             color="black", fill="darkorchid3", shape=22, size=2.5)+
  geom_point(aes(y=3.5,x=9.4),
             color="black", fill="deepskyblue3", shape=23, size=2.5)+
  geom_point(aes(y=2.6,x=9.7),
             color="black", fill="gold2", shape=24, size=2.5)+
  #matched points treatment
  geom_point(aes(y=3.8,x=10.95),
             color="black", fill="deepskyblue3", shape=22, size=2.5)+
  geom_point(aes(y=2.8,x=10.7),
             color="black", fill="darkorchid3", shape=23, size=2.5)+
  geom_point(aes(y=2.7,x=10.3),
             color="black", fill="gold2", shape=24, size=2.5)



### 8

ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="light grey")+geom_point(fill="light grey", color="light grey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_vline(aes(xintercept=9),linetype='dashed', colour="dark grey", size=1) +
  geom_vline(aes(xintercept=11),linetype='dashed', colour="dark grey", size=1) + 
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "BW1", hjust=0.5, size=4, col="black")+
  geom_segment(aes(x = 9.4, y = 11, xend = 9, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 11, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="darkorchid3", shape=22, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="deepskyblue3", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="gold2", shape=24, size=3)+
  #Average
  geom_segment(aes(x = 9.5, y = 3.07, xend = 10, yend = 3.07), col="black", size=0.8,lty=6) +
  geom_segment(aes(x = 10, y = 3.1, xend = 10.5, yend = 3.1), col="black", size=0.8,lty=6) +
  geom_segment(aes(x = 10, y = 3.07, xend = 10, yend = 3.1), col="#eb811b", size=1) +
  geom_point(aes(y=3.07,x=9.5),
             color="red", fill="red", size=3)+
  geom_point(aes(y=3.1,x=10.5),
             color="blue", fill="blue", size=3)+
  annotate("rect", xmin = 16, xmax = 20.2, ymin = -1.2, ymax = 0.5,
           color="black", fill="white") +
  annotate("text", x = 16.8, y = 0, label = "Avg. matched T", hjust=0, vjust=0.5, size=5)+
  annotate("text", x = 16.8, y = -0.7, label = "Avg. matched C", hjust=0, vjust=0.5, size=5)+
  geom_point(aes(y=0,x=16.3),
             color="red", fill="red", size=3)+
  geom_point(aes(y=-0.7,x=16.3),
             color="blue", fill="blue", size=3)

### 9

ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_vline(aes(xintercept=9),linetype='dashed', colour="light grey", size=1) +
  geom_vline(aes(xintercept=11),linetype='dashed', colour="light grey", size=1) + 
  
  geom_vline(aes(xintercept=8),linetype='dashed', colour="dark grey", size=1) +
  geom_vline(aes(xintercept=12),linetype='dashed', colour="dark grey", size=1) + 
  
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "BW2", hjust=0.5, size=4, col="black")+
  geom_segment(aes(x = 9.4, y = 11, xend = 8, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 12, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="yellow2", shape=25, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="slateblue3", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="turquoise3", shape=21, size=3)



### 10

ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="light grey")+geom_point(fill="light grey", color="light grey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_vline(aes(xintercept=9),linetype='dashed', colour="light grey", size=1) +
  geom_vline(aes(xintercept=11),linetype='dashed', colour="light grey", size=1) + 
  
  geom_vline(aes(xintercept=8),linetype='dashed', colour="dark grey", size=1) +
  geom_vline(aes(xintercept=12),linetype='dashed', colour="dark grey", size=1) + 
  
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "BW2", hjust=0.5, size=4, col="black")+
  geom_segment(aes(x = 9.4, y = 11, xend = 8, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 12, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="yellow2", shape=25, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="slateblue3", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="turquoise3", shape=21, size=3)+
  #Matched obs (control)
  geom_point(aes(y=3.5,x=8.2),
             color="black", fill="yellow2", shape=25, size=2.5)+
  geom_point(aes(y=4.05,x=8.9),
             color="black", fill="slateblue3", shape=23, size=2.5)+
  geom_point(aes(y=2.5,x=8.4),
             color="black", fill="turquoise3", shape=21, size=2.5) +
  #matchd obs (treatment)
  geom_point(aes(y=3.2,x=11.1),
             color="black", fill="yellow2", shape=25, size=2.5)+
  geom_point(aes(y=4.5,x=11.9),
             color="black", fill="slateblue3", shape=23, size=2.5)+
  geom_point(aes(y=2.6,x=11.5),
             color="black", fill="turquoise3", shape=21, size=2.5)


### 11

ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="light grey")+geom_point(fill="light grey", color="light grey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_vline(aes(xintercept=9),linetype='dashed', colour="light grey", size=1) +
  geom_vline(aes(xintercept=11),linetype='dashed', colour="light grey", size=1) + 
  
  geom_vline(aes(xintercept=8),linetype='dashed', colour="dark grey", size=1) +
  geom_vline(aes(xintercept=12),linetype='dashed', colour="dark grey", size=1) + 
  
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "BW2", hjust=0.5, size=4, col="black")+
  geom_segment(aes(x = 9.4, y = 11, xend = 8, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 12, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="yellow2", shape=25, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="slateblue3", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="turquoise3", shape=21, size=3)+
  #Average
  geom_segment(aes(x = 9, y = 3.35, xend = 10, yend = 3.35), col="black", size=0.8,lty=6) +
  geom_segment(aes(x = 10, y = 3.4, xend = 11, yend = 3.4), col="black", size=0.8,lty=6) +
  geom_segment(aes(x = 10, y = 3.35, xend = 10, yend = 3.4), col="#eb811b", size=1) +
  geom_point(aes(y=3.35,x=9),
             color="red", fill="red", size=3)+
  geom_point(aes(y=3.4,x=11),
             color="blue", fill="blue", size=3)+
  annotate("rect", xmin = 16, xmax = 20.2, ymin = -1.2, ymax = 0.5,
           color="black", fill="white") +
  annotate("text", x = 16.8, y = 0, label = "Avg. matched T", hjust=0, vjust=0.5, size=5)+
  annotate("text", x = 16.8, y = -0.7, label = "Avg. matched C", hjust=0, vjust=0.5, size=5)+
  geom_point(aes(y=0,x=16.3),
             color="red", fill="red", size=3)+
  geom_point(aes(y=-0.7,x=16.3),
             color="blue", fill="blue", size=3)




### 12

ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="light grey")+geom_point(fill="light grey", color="light grey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_vline(aes(xintercept=9),linetype='dashed', colour="light grey", size=1) +
  geom_vline(aes(xintercept=11),linetype='dashed', colour="light grey", size=1) + 
  
  geom_vline(aes(xintercept=8),linetype='dashed', colour="light grey", size=1) +
  geom_vline(aes(xintercept=12),linetype='dashed', colour="light grey", size=1) + 
  
  geom_vline(aes(xintercept=7),linetype='dashed', colour="dark grey", size=1) +
  geom_vline(aes(xintercept=13),linetype='dashed', colour="dark grey", size=1) + 
  
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "BW3", hjust=0.5, size=4, col="black")+
  geom_segment(aes(x = 9.4, y = 11, xend = 7, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 13, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="deeppink2", shape=22, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="mediumseagreen", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="sienna2", shape=24, size=3)+

  #Average
  geom_segment(aes(x = 8.5, y = 3.2, xend = 10, yend = 3.2), col="black", size=0.8,lty=6) +
  geom_segment(aes(x = 10, y = 3.6, xend = 11.5, yend = 3.6), col="black", size=0.8,lty=6) +
  geom_segment(aes(x = 10, y = 3.2, xend = 10, yend = 3.6), col="#eb811b", size=1) +
  geom_point(aes(y=3.2,x=8.5),
             color="red", fill="red", size=3)+
  geom_point(aes(y=3.6,x=11.5),
             color="blue", fill="blue", size=3)+
  annotate("rect", xmin = 16, xmax = 20.2, ymin = -1.2, ymax = 0.5,
           color="black", fill="white") +
  annotate("text", x = 16.8, y = 0, label = "Avg. matched T", hjust=0, vjust=0.5, size=5)+
  annotate("text", x = 16.8, y = -0.7, label = "Avg. matched C", hjust=0, vjust=0.5, size=5)+
  geom_point(aes(y=0,x=16.3),
             color="red", fill="red", size=3)+
  geom_point(aes(y=-0.7,x=16.3),
             color="blue", fill="blue", size=3)


### 13

ggplot(df_pre,aes(y=Y,x=xaxisTime, label=NA), fill="light grey")+geom_point(fill="light grey", color="light grey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  geom_vline(aes(xintercept=9),linetype='dashed', colour="light grey", size=1) +
  geom_vline(aes(xintercept=11),linetype='dashed', colour="light grey", size=1) + 
  
  geom_vline(aes(xintercept=8),linetype='dashed', colour="deeppink2", size=1.5) +
  geom_vline(aes(xintercept=12),linetype='dashed', colour="deeppink2", size=1.5) + 
  
  geom_vline(aes(xintercept=7),linetype='dashed', colour="light grey", size=1) +
  geom_vline(aes(xintercept=13),linetype='dashed', colour="light grey", size=1) + 
  
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "H*", hjust=0.5, size=5, col="black", fontface=2)+
  geom_segment(aes(x = 9.4, y = 11, xend = 8, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 12, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="yellow2", shape=25, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="slateblue3", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="turquoise3", shape=21, size=3)



### 14

ggplot(df_post,aes(y=Y,x=xaxisTime, label=NA), fill="slategrey")+geom_point(fill="slategrey", color="slategrey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+

  geom_vline(aes(xintercept=8),linetype='dashed', colour="deeppink2", size=1.5) +
  geom_vline(aes(xintercept=12),linetype='dashed', colour="deeppink2", size=1.5) + 
  
  xlab("Running variable") + ylab("Outcome")+ggtitle("Post-intervention")+
  
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "H*", hjust=0.5, size=5, col="black", fontface=2)+
  geom_segment(aes(x = 9.4, y = 11, xend = 8, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 12, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="yellow2", shape=25, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="slateblue3", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="turquoise3", shape=21, size=3)


### 15

ggplot(df_post,aes(y=Y,x=xaxisTime, label=NA), fill="light grey")+geom_point(fill="light grey", color="light grey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  
  geom_vline(aes(xintercept=8),linetype='dashed', colour="deeppink2", size=1.5) +
  geom_vline(aes(xintercept=12),linetype='dashed', colour="deeppink2", size=1.5) + 
  
  xlab("Running variable") + ylab("Outcome")+ggtitle("Post-intervention")+
  
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "H*", hjust=0.5, size=5, col="black", fontface=2)+
  geom_segment(aes(x = 9.4, y = 11, xend = 8, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 12, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="yellow2", shape=25, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="slateblue3", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="turquoise3", shape=21, size=3) +
  
  #Matched obs (control)
  geom_point(aes(y=3.8,x=8.1),
             color="black", fill="yellow2", shape=25, size=2.5)+
  geom_point(aes(y=4.05,x=8.9),
             color="black", fill="slateblue3", shape=23, size=2.5)+
  geom_point(aes(y=3.7,x=9.9),
             color="black", fill="turquoise3", shape=21, size=2.5) +
  #matchd obs (treatment)
  geom_point(aes(y=5.8,x=10.7),
             color="black", fill="yellow2", shape=25, size=2.5)+
  geom_point(aes(y=6.5,x=10.05),
             color="black", fill="slateblue3", shape=23, size=2.5)+
  geom_point(aes(y=5.9,x=11.9),
             color="black", fill="turquoise3", shape=21, size=2.5)


### 16

ggplot(df_post,aes(y=Y,x=xaxisTime, label=NA), fill="light grey")+geom_point(fill="light grey", color="light grey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  
  geom_vline(aes(xintercept=8),linetype='dashed', colour="deeppink2", size=1.5) +
  geom_vline(aes(xintercept=12),linetype='dashed', colour="deeppink2", size=1.5) + 
  
  xlab("Running variable") + ylab("Outcome")+ggtitle("Post-intervention")+
  
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "H*", hjust=0.5, size=5, col="black", fontface=2)+
  geom_segment(aes(x = 9.4, y = 11, xend = 8, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 12, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="yellow2", shape=25, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="slateblue3", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="turquoise3", shape=21, size=3) +
  
  #Average
  geom_segment(aes(x = 9, y = 3.85, xend = 10, yend = 3.85), col="black", size=0.8,lty=6) +
  geom_segment(aes(x = 10, y = 6.07, xend = 11, yend = 6.07), col="black", size=0.8,lty=6) +
  geom_segment(aes(x = 10, y = 3.85, xend = 10, yend = 6.07), col="#eb811b", size=1) +
  geom_point(aes(y=3.85,x=9),
             color="red", fill="red", size=3)+
  geom_point(aes(y=6.07,x=11),
             color="blue", fill="blue", size=3)+
  annotate("rect", xmin = 16, xmax = 20.2, ymin = -1.2, ymax = 0.5,
           color="black", fill="white") +
  annotate("text", x = 16.8, y = 0, label = "Avg. matched T", hjust=0, vjust=0.5, size=5)+
  annotate("text", x = 16.8, y = -0.7, label = "Avg. matched C", hjust=0, vjust=0.5, size=5)+
  geom_point(aes(y=0,x=16.3),
             color="red", fill="red", size=3)+
  geom_point(aes(y=-0.7,x=16.3),
             color="blue", fill="blue", size=3)



### 17

ggplot(df_post,aes(y=Y,x=xaxisTime, label=NA), fill="light grey")+geom_point(fill="light grey", color="light grey")+
  geom_vline(aes(xintercept=10),linetype='dashed') + theme_bw()+
  
  geom_vline(aes(xintercept=8),linetype='dashed', colour="deeppink2", size=1.5) +
  geom_vline(aes(xintercept=12),linetype='dashed', colour="deeppink2", size=1.5) + 
  
  xlab("Running variable") + ylab("Outcome")+ggtitle("Post-intervention")+
  
  annotate("rect", xmin = 9.1, xmax = 10.9, ymin = 10.6, ymax = 11.4,
           color="white", fill="white") +
  annotate("text", x = 10, y = 11, label = "H*", hjust=0.5, size=5, col="black", fontface=2)+
  geom_segment(aes(x = 9.4, y = 11, xend = 8, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  geom_segment(aes(x = 10.6, y = 11, xend = 12, yend = 11), col="black", size=0.7,
               arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "open"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        title=element_text(size=15,hjust=0.5, face="bold"))+
  annotate("rect", xmin = 0, xmax = 3, ymin = 10, ymax = 11.5,
           color="black", fill="white") +
  annotate("text", x = 1.5, y = 11, label = "Template", hjust=0.5, vjust=0, size=5)+
  geom_point(aes(y=10.5,x=0.5),
             color="black", fill="yellow2", shape=25, size=3)+
  geom_point(aes(y=10.5,x=1.5),
             color="black", fill="slateblue3", shape=23, size=3)+
  geom_point(aes(y=10.5,x=2.5),
             color="black", fill="turquoise3", shape=21, size=3) +
  
  #Average
  geom_segment(aes(x = 9, y = 3.85, xend = 10, yend = 3.85), col="black", size=0.8,lty=6) +
  geom_segment(aes(x = 10, y = 6.07, xend = 11, yend = 6.07), col="black", size=0.8,lty=6) +
  geom_segment(aes(x = 10, y = 3.85, xend = 10, yend = 6.07), col="#eb811b", size=1) +
  geom_point(aes(y=3.85,x=9),
             color="red", fill="red", size=3)+
  geom_point(aes(y=6.07,x=11),
             color="blue", fill="blue", size=3)+
  annotate("text", x = 8.5, y = 5, label = "TOT", hjust=0, size=5, col="#eb811b", fontface=2)+
  annotate("rect", xmin = 16, xmax = 20.2, ymin = -1.2, ymax = 0.5,
           color="black", fill="white") +
  annotate("text", x = 16.8, y = 0, label = "Avg. matched T", hjust=0, vjust=0.5, size=5)+
  annotate("text", x = 16.8, y = -0.7, label = "Avg. matched C", hjust=0, vjust=0.5, size=5)+
  geom_point(aes(y=0,x=16.3),
             color="red", fill="red", size=3)+
  geom_point(aes(y=-0.7,x=16.3),
             color="blue", fill="blue", size=3)
