## Self contained code

# Title: Simulations
# Paper: How far is too far?
# Authors: M. Bennett
# Date: 08-06-19

## Load libraries
library(designmatch)
library(Rglpk)
library(gurobi)
library(xtable)
library(dplyr)
library(readstata13)
library(rdrobust)
library(exact2x2)
library(simstudy)
library(exact2x2)
library(foreign)
library(rbounds)
library(stats)
library(sensitivitymv)
library(DTComPair)
library(Rlab)
library(nprobust)
library(KernSmooth)
library(lpdensity)
library(Smisc)
library(ggplot2)
library(firasans)

##############################################################################################

#sim 3, 33 corr low sample large

githubURL <- "https://github.com/maibennett/presentations/raw/main/content/presentations/RD/IC_20210312/data/example_grd.RData"
load(url(githubURL))

col_text = "#23373B"

#palette v2
col1 = "#8507C4"
#col2 = "#28A4FF"
#col3 = "#FFDB00"

#palette v3
#col1 = "#8A1C7C"
#col1 = "#8A1C8D"
col2 = "#EE6352"
col3 = "#FFCF4F"

par(mfrow=c(1,1),mar=c(4,3.5,1.5,2))

set.seed(1)

random = sample(1:nrow(data),2000,replace=FALSE)

data2 = data[random,]


### Plot 1

ggplot(data2[data2[,tvar]==0,],aes(y=y,x=running_var, label=NA), fill="white")+geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=1) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1.5) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1.5)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1.5)+
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  #  annotate("pointrage",x=-600,y=41,xend=-600,yend=41,shape=21,size=2,col="dark grey") +
  #  annotate("text",x=13.05,y=8.4,label = "H*", size=5.5, colour=col_text)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-60,60) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=35,ymin=35,ymax=35,shape=21,size=1,col="slategrey",fill="white",alpha=1)+
  annotate("text",x=-550,y=35,label = "Observations", size=5, colour=col_text,hjust = 0,family=font_fsc)




## Plot 2

#RUN THE FIRST ITERATION

#d_match = grd$d_rematch

random_t = sample(1:nrow(d_match),600,replace=FALSE)
d_match2 = d_match[random_t,]


ggplot(data2[data2[,tvar]==0,],aes(y=y,x=running_var, label=NA), fill="white")+
  geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2[d_match2[,running_var]>=bandwidth[1] & bandwidth[2]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color='light grey', size=2, shape=21)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=1) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1.5) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1.5)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1.5)+
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  #  annotate("pointrage",x=-600,y=41,xend=-600,yend=41,shape=21,size=2,col="dark grey") +
  #  annotate("text",x=13.05,y=8.4,label = "H*", size=5.5, colour=col_text)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-60,60) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=34,ymin=34,ymax=34,shape=21,size=1,col=col2,fill=col2,alpha=1)+
  annotate("text",x=-550,y=34,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)



ggplot(data2[data2[,tvar]==0,],aes(y=y,x=running_var, label=NA), fill="white")+
  geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2[d_match2[,running_var]>=bandwidth[1] & bandwidth[2]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_point(data=d_match2[d_match2[,running_var]>=grid[4] & grid[5]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color='light grey', size=2, shape=21)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=1) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1.5) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1.5)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1.5)+
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  #  annotate("pointrage",x=-600,y=41,xend=-600,yend=41,shape=21,size=2,col="dark grey") +
  #  annotate("text",x=13.05,y=8.4,label = "H*", size=5.5, colour=col_text)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-60,60) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=34,ymin=34,ymax=34,shape=21,size=1,col=col2,fill=col2,alpha=0.2)+
  annotate("text",x=-550,y=34,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)




ggplot(data2[data2[,tvar]==0,],aes(y=y,x=running_var, label=NA), fill="white")+
  geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2[d_match2[,running_var]>=bandwidth[1] & bandwidth[2]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_point(data=d_match2[d_match2[,running_var]>=grid[4] & grid[5]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_point(data=d_match2[d_match2[,running_var]>=grid[3] & grid[4]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color='light grey', size=2, shape=21,alpha=1)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=1) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1.5) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1.5)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1.5)+
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  #  annotate("pointrage",x=-600,y=41,xend=-600,yend=41,shape=21,size=2,col="dark grey") +
  #  annotate("text",x=13.05,y=8.4,label = "H*", size=5.5, colour=col_text)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-60,60) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=34,ymin=34,ymax=34,shape=21,size=1,col=col2,fill=col2,alpha=0.2)+
  annotate("text",x=-550,y=34,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)



ggplot(data2[data2[,tvar]==0,],aes(y=y,x=running_var, label=NA), fill="white")+
  geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2[d_match2[,running_var]>=bandwidth[1] & bandwidth[2]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_point(data=d_match2[d_match2[,running_var]>=grid[4] & grid[5]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_point(data=d_match2[d_match2[,running_var]>=grid[3] & grid[4]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_point(data=d_match2[d_match2[,running_var]>=grid[2] & grid[3]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color='light grey', size=2, shape=21,alpha=1)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=1) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1.5) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1.5)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1.5)+
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  #  annotate("pointrage",x=-600,y=41,xend=-600,yend=41,shape=21,size=2,col="dark grey") +
  #  annotate("text",x=13.05,y=8.4,label = "H*", size=5.5, colour=col_text)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-60,60) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=34,ymin=34,ymax=34,shape=21,size=1,col=col2,fill=col2,alpha=0.2)+
  annotate("text",x=-550,y=34,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)




ggplot(d_match2,aes(y=y,x=running_var, label=NA), fill=col2)+
#  geom_rect(aes(xmin = grid[2], xmax = grid[9], ymin = -Inf, ymax = Inf),
#            fill = col2, alpha = 0.03) +
  geom_ribbon(aes(x=running_var,ymin = -60, ymax = 60),fill = col2, alpha = 0.3) +
  geom_point(fill=col2, color='light grey', size=2, shape=21)+
  geom_point(data=data2[data2[,tvar]==0,],aes(y=y,x=running_var, label=NA), fill="white",col="slategrey",shape=21,alpha=0.5)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=0.6) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1.1) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1)+
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-60,60) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=35,ymin=35,ymax=35,shape=21,size=1,col=col2,fill=col2,alpha=1)+
  annotate("text",x=-550,y=35,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)+
  annotate("rect",xmin=-600,xmax=-580,ymin=23,ymax=27,col="dark gray",fill=col2,alpha=0.3)+
  annotate("text",x=-550,y=25,label = "Overlap region", size=5, colour=col_text,hjust = 0,family=font_fsc)
  


### Plot 4

data_eval = as.data.frame(cbind(eval,est,cil,ciu))
names(data_eval) = c("running_var","est","cil","ciu")

data3 = data2[data2[,tvar]==0,]

library(gtools)
data3 = smartbind(data3,data_eval)

ggplot(data3,aes(y=y,x=running_var, label=NA), fill="white")+
  #geom_rect(aes(xmin = grid[5], xmax = grid[6], ymin = -Inf, ymax = Inf),
  #          fill = "light grey", alpha = 0.03) +
  geom_ribbon(data=data3[data3$running_var>=grid[5] & data3$running_var<=grid[6],],
              aes(ymin=-30,ymax=30),fill="light grey",alpha=0.5) +
  #geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2,
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.3)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=0.6) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1)+
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=est),col=col_text,lwd=1) +
  geom_ribbon(data=data3[!is.na(data3$est),],aes(x=running_var,ymin=cil,ymax=ciu),fill=col1,alpha=0.3) +
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=cil), col=col_text, lty=4) +
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=ciu), col=col_text, lty=4) +
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-30,30) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=18,ymin=18,ymax=18,shape=21,size=1,col=col2,fill=col2,alpha=0.2)+
  annotate("text",x=-550,y=18,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)+
  annotate("rect",xmin=-600,xmax=-580,ymin=12.3,ymax=13.9,col="dark gray",fill=col1,alpha=0.3)+
  annotate("text",x=-550,y=13.4,label = "Loc. Poly. 95% CI", size=5, colour=col_text,hjust = 0,family=font_fsc) +
  annotate("rect",xmin=-600,xmax=-580,ymin=8,ymax=9.5,col="dark gray",fill="light grey",alpha=0.5)+
  annotate("text",x=-550,y=9,label = "Template 1 region", size=5, colour=col_text,hjust = 0,family=font_fsc)


ggplot(data3,aes(y=y,x=running_var, label=NA), fill="white")+
  #geom_rect(aes(xmin = grid[5], xmax = grid[6], ymin = -Inf, ymax = Inf),
  #          fill = "light grey", alpha = 0.03) +
  geom_ribbon(data=data3[data3$running_var>=grid[5] & data3$running_var<=grid[6],],
              aes(ymin=-30,ymax=30),fill="light grey",alpha=0.5) +
  #geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2,
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.3)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=0.6) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=optimal_bandwidth[1]),linetype='dashed',col=col1,lwd=1.1) + 
  geom_vline(aes(xintercept=optimal_bandwidth[2]),linetype='dashed',col=col1,lwd=1.1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col=col1), size=1)+
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=est),col=col_text,lwd=1) +
  geom_ribbon(data=data3[!is.na(data3$est),],aes(x=running_var,ymin=cil,ymax=ciu),fill=col1,alpha=0.3) +
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=cil), col=col_text, lty=4) +
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=ciu), col=col_text, lty=4) +
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black',col1),
                       labels = c("Grid","Cutoff","H1"),
                       guide = "legend")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-30,30) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=15,ymin=15,ymax=15,shape=21,size=1,col=col2,fill=col2,alpha=0.2)+
  annotate("text",x=-550,y=15,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)+
  annotate("rect",xmin=-600,xmax=-580,ymin=10,ymax=11.5,col="dark gray",fill=col1,alpha=0.3)+
  annotate("text",x=-550,y=10.5,label = "Loc. Poly. 95% CI", size=5, colour=col_text,hjust = 0,family=font_fsc) +
  annotate("rect",xmin=-600,xmax=-580,ymin=5.2,ymax=6.8,col="dark gray",fill="light grey",alpha=0.5)+
  annotate("text",x=-550,y=6,label = "Template region", size=5, colour=col_text,hjust = 0,family=font_fsc)



## RUN ALL ITERATIONS
githubURL <- "https://github.com/maibennett/presentations/raw/main/content/presentations/RD/IC_20210312/data/example_grd_all.RData"
load(url(githubURL))

random_t = sample(1:nrow(d_match),600,replace=FALSE)
d_match2 = d_match[random_t,]

data_eval_1 = data_eval

data_eval = as.data.frame(cbind(eval,est,cil,ciu))
names(data_eval) = c("running_var","est","cil","ciu")

data_eval_1_trim = data_eval_1[data_eval_1$running_var<max(data_eval$running_var),]
  
data3 = data2[data2[,tvar]==0,]

library(gtools)
#data3 = smartbind(data3,data_eval)
data3 = smartbind(data3,data_eval_1_trim)



ggplot(data3,aes(y=y,x=running_var, label=NA), fill="white")+
  #geom_rect(aes(xmin = bandwidth_opt[1]+15, xmax = bandwidth_opt[2], ymin = -Inf, ymax = Inf),
  #          fill = col3, alpha = 0.003) +
  geom_ribbon(data=data3[data3$running_var>=bandwidth_opt[1]+15 & data3$running_var<=grid[6],],
              aes(ymin=-30,ymax=30),fill="light grey",alpha=0.5) +
 # geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2,
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=0.6) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=optimal_bandwidth[1]+15),linetype='dashed',col=col1,lwd=1.1) + 
  geom_vline(aes(xintercept=optimal_bandwidth[2]),linetype='dashed',col=col1,lwd=1.1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col=col1), size=1.5)+
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=est),col=col_text,lwd=1) +
  geom_ribbon(data=data3[!is.na(data3$est),],aes(x=running_var,ymin=cil,ymax=ciu),fill=col1,alpha=0.3) +
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=cil), col=col_text, lty=4) +
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=ciu), col=col_text, lty=4) +
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black',col1),
                       labels = c("Grid","Cutoff","H*"),
                       guide = "legend")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-30,30) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=15,ymin=15,ymax=15,shape=21,size=1,col=col2,fill=col2,alpha=0.2)+
  annotate("text",x=-550,y=15,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)+
  annotate("rect",xmin=-600,xmax=-580,ymin=10,ymax=11.5,col="dark gray",fill=col1,alpha=0.3)+
  annotate("text",x=-550,y=10.5,label = "Loc. Poly. 95% CI", size=5, colour=col_text,hjust = 0,family=font_fsc) +
  annotate("rect",xmin=-600,xmax=-580,ymin=5.2,ymax=6.8,col="dark gray",fill="light grey",alpha=0.5)+
  annotate("text",x=-550,y=6,label = "Template region", size=5, colour=col_text,hjust = 0,family=font_fsc)

