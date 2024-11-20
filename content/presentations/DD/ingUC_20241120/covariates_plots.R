#Clear memory
rm(list = ls())

#Clear the console
cat("\014")

library(ggplot2)
library(dplyr)
library(mvtnorm)
library(tidyverse)
library(designmatch)
library(hrbrthemes)
library(ggridges)
library(multiwayvcov)
library(lmtest)
library(sandwich)
library(modelsummary)
library(MatchIt)

#Turn off scientific notation (turn back on with 0)
set.seed(100)

### Simulations

col1 = "#5601A4"
col2 = "#FCCE25"
  
s = 1

t_out = c(2,4,6,8,10)  

########### Scenario 1: Time-invariant covariate effect
i = 1
set.seed(as.numeric(paste0(s,i)))

n = 1000
max.time = 10
trt.time = 6

d <- expand.grid(id = 1:n, t = 1:max.time) %>%
  arrange(id,t) %>% group_by(id) %>%
  mutate(u=rnorm(1,0,sd=0.25), # random intercept
         p.trt=0.5, # probability of treatment
         a=rbinom(1, 1, p.trt), # treatment
         x=rnorm(1, mean = 1.5 - 0.5*a, sd = 1.5 - 0.5*a),
         p=I(t >= trt.time), # indicator of post-treatment period
         treated=I(p == 1 & a == 1), # treated indicator
  ) %>%
  ungroup()

d <- d %>% mutate(err=rnorm(n*max.time),
                  y = 1 + x + a + u + treated +
                    ((t - 2.5)^2)/10 + err) %>%
  group_by(id) %>%
  mutate(y.diff = y - lag(y)) %>%
  ungroup()


d = d[!(d$t %in% t_out),]

g1 = ggplot(data = d, aes(x = x, y = factor(t), fill = factor(a), color=factor(a)),fill="white") +
  stat_density_ridges(alpha=0.3, scale = .9,lwd=1.2, quantile_lines = TRUE, quantiles = 0.5) + 
  scale_fill_manual(values=c("0"=col1, "1"=col2), labels = c("C","T")) +
  scale_color_manual(values=c("0"=col1, "1"=col2), labels = c("C","T")) +
  theme_bw()+
  theme_ipsum_rc(plot_title_face = "bold") + #plain 
  xlab("") + ylab("Time") + ggtitle("Scenario 1")+
  
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(axis.title.x = element_text(size=14),#margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=18),#margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=10),legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.background = element_rect(fill="white",colour ="white"),
        title = element_text(size=9))


########### Scenario 2: Time-varying covariate effect

i = 2
set.seed(as.numeric(paste0(s,i)))

n = 1000
max.time = 10
trt.time = 6

d <- expand.grid(id = 1:n, t = 1:max.time) %>%
  arrange(id,t) %>% group_by(id) %>%
  mutate(u=rnorm(1,0,sd=0.25), # random intercept
         p.trt=0.5, # probability of treatment
         a=rbinom(1, 1, p.trt), # treatment
         x=rnorm(1, mean = 1.5 - 0.5*a, sd = 1.5 - 0.5*a),
         p=I(t >= trt.time), # indicator of post-treatment period
         treated=I(p == 1 & a == 1), # treated indicator
  ) %>%
  ungroup()

d <- d %>% mutate(err=rnorm(n*max.time),
                  y = 1 + x + a + u + treated +
                    ((t - 2.5)^2)/10 + x*t/10 + err) %>%
  group_by(id) %>%
  mutate(y.diff = y - lag(y)) %>%
  ungroup()

d = d[!(d$t %in% t_out),]

g2 = ggplot(data = d, aes(x = x, y = factor(t), fill = factor(a), color=factor(a)),fill="white") +
  stat_density_ridges(alpha=0.3, scale = .9,lwd=1.2, quantile_lines = TRUE, quantiles = 0.5) + 
  scale_fill_manual(values=c("0"=col1, "1"=col2), labels = c("C","T")) +
  scale_color_manual(values=c("0"=col1, "1"=col2), labels = c("C","T")) +
  theme_bw()+
  theme_ipsum_rc(plot_title_face = "bold") + #plain 
  xlab("") + ylab("") + ggtitle("Scenario 2")+
  
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(axis.title.x = element_text(size=14),#margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=14),#margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=10),legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.background = element_rect(fill="white",colour ="white"),
        title = element_text(size=14))


########### Scenario 3: Treatment independent covariate

i = 3
set.seed(as.numeric(paste0(s,i)))

n = 1000
max.time = 10
trt.time = 6

d <- expand.grid(id = 1:n, t = 1:max.time) %>%
  arrange(id,t) %>% group_by(id) %>%
  mutate(u=rnorm(1,0,sd=0.25), # random intercept
         p.trt=0.5, # probability of treatment
         a=rbinom(1, 1, p.trt), # treatment
         x=rnorm(1, mean = 1, sd = 1),
         p=I(t >= trt.time), # indicator of post-treatment period
         treated=I(p == 1 & a == 1), # treated indicator
  ) %>%
  ungroup()

d <- d %>% mutate(err=rnorm(n*max.time),
                  y = 1 + x + a + u + treated +
                    ((t - 2.5)^2)/10 + x*t/10 + err) %>%
  group_by(id) %>%
  mutate(y.diff = y - lag(y)) %>%
  ungroup()

d = d[!(d$t %in% t_out),]

g3 = ggplot(data = d, aes(x = x, y = factor(t), fill = factor(a), color=factor(a)),fill="white") +
  stat_density_ridges(alpha=0.3, scale = .9,lwd=1.2, quantile_lines = TRUE, quantiles = 0.5) + 
  scale_fill_manual(values=c("0"=col1, "1"=col2), labels = c("C","T")) +
  scale_color_manual(values=c("0"=col1, "1"=col2), labels = c("C","T")) +
  theme_bw()+
  theme_ipsum_rc(plot_title_face = "bold") + #plain 
  xlab("X") + ylab("") + ggtitle("Scenario 3")+
  
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(axis.title.x = element_text(size=18),#margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=14),#margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=10),legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.background = element_rect(fill="white",colour ="white"),
        title = element_text(size=14))


########### Scenario 4: Parallel evolution

i = 4
set.seed(as.numeric(paste0(s,i)))

n = 1000
max.time = 10
trt.time = 6

d <- expand.grid(id = 1:n, t = 1:max.time) %>%
  arrange(id,t) %>% group_by(id) %>%
  mutate(u=rnorm(1,0,sd=0.25), # random intercept
         p.trt=0.5, # probability of treatment
         a=rbinom(1, 1, p.trt), # treatment
         x=rnorm(1, mean = 1.5 - 0.5*a, sd = 1.5 - 0.5*a),
         p=I(t >= trt.time), # indicator of post-treatment period
         treated=I(p == 1 & a == 1), # treated indicator
         x=ifelse(t>=2, lag(x, 1) + (t-1)/4 *
                    rnorm(1, mean = 1, sd = 0.1), x)
  ) %>%
  ungroup()

d <- d %>% mutate(err=rnorm(n*max.time),
                  y = 1 + x + a + u + treated +
                    ((t - 2.5)^2)/10 + x*t/10 + err) %>%
  group_by(id) %>%
  mutate(y.diff = y - lag(y)) %>%
  ungroup()

d = d[!(d$t %in% t_out),]

g4 = ggplot(data = d, aes(x = x, y = factor(t), fill = factor(a), color=factor(a)),fill="white") +
  stat_density_ridges(alpha=0.3, scale = .9,lwd=1.2, quantile_lines = TRUE, quantiles = 0.5) + 
  scale_fill_manual(values=c("0"=col1, "1"=col2), labels = c("C","T")) +
  scale_color_manual(values=c("0"=col1, "1"=col2), labels = c("C","T")) +
  theme_bw()+
  theme_ipsum_rc(plot_title_face = "bold") + #plain 
  xlab("") + ylab("Time") + ggtitle("Scenario 4")+
  
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(axis.title.x = element_text(size=14),#margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=18),#margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=10),legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.background = element_rect(fill="white",colour ="white"),
        title = element_text(size=14))


########### Scenario 5: Evolution differs by group

i = 5
set.seed(as.numeric(paste0(s,i)))

n = 1000
max.time = 10
trt.time = 6

d <- expand.grid(id = 1:n, t = 1:max.time) %>%
  arrange(id,t) %>% group_by(id) %>%
  mutate(u=rnorm(1,0,sd=0.25), # random intercept
         p.trt=0.5, # probability of treatment
         a=rbinom(1, 1, p.trt), # treatment
         x=rnorm(1, mean = 1.5 - 0.5*a, sd = 1.5 - 0.5*a),
         p=I(t >= trt.time), # indicator of post-treatment period
         treated=I(p == 1 & a == 1), # treated indicator
         x=ifelse(t>=2 & a==1, lag(x, 1) + (t-1)/10 *
                    rnorm(1, mean = 1, sd = 0.1), x),
         x=ifelse(t>=2 & a==0, lag(x, 1) - (t-1)/10 *
                    rnorm(1, mean = 1, sd = 0.1), x)
  ) %>%
  ungroup()

d <- d %>% mutate(err=rnorm(n*max.time),
                  y = 1 + x + a + u + treated +
                    ((t - 2.5)^2)/10 + x*t/10 + err) %>%
  group_by(id) %>%
  mutate(y.diff = y - lag(y)) %>%
  ungroup()

d = d[!(d$t %in% t_out),]

g5 = ggplot(data = d, aes(x = x, y = factor(t), fill = factor(a), color=factor(a)),fill="white") +
  stat_density_ridges(alpha=0.3, scale = .9,lwd=1.2, quantile_lines = TRUE, quantiles = 0.5) + 
  scale_fill_manual(values=c("0"=col1, "1"=col2), labels = c("C","T")) +
  scale_color_manual(values=c("0"=col1, "1"=col2), labels = c("C","T")) +
  theme_bw()+
  theme_ipsum_rc(plot_title_face = "bold") + #plain 
  xlab("") + ylab("") + ggtitle("Scenario 5")+
  
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(axis.title.x = element_text(size=14),#margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=14),#margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=10),legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.background = element_rect(fill="white",colour ="white"),
        title = element_text(size=14))



########### Scenario 6: Evolution diverges in post

i = 6
set.seed(as.numeric(paste0(s,i)))

n = 1000
max.time = 10
trt.time = 6

d <- expand.grid(id = 1:n, t = 1:max.time) %>%
  arrange(id,t) %>% group_by(id) %>%
  mutate(u=rnorm(1,0,sd=0.25), # random intercept
         p.trt=0.5, # probability of treatment
         a=rbinom(1, 1, p.trt), # treatment
         x=rnorm(1, mean = 1.5 - 0.5*a, sd = 1.5 - 0.5*a),
         p=I(t >= trt.time), # indicator of post-treatment period
         treated=I(p == 1 & a == 1), # treated indicator
         x=ifelse(t>=2, lag(x, 1) + (t-1)/10 *
                    rnorm(1, mean = 1, sd = 0.1) -
                    a*p*t/20, x)
  ) %>%
  ungroup()

d <- d %>% mutate(err=rnorm(n*max.time),
                  y = 1 + x + a + u + treated +
                    ((t - 2.5)^2)/10 + x*t/10 + err) %>%
  group_by(id) %>%
  mutate(y.diff = y - lag(y)) %>%
  ungroup()


effect = mean(d$y[d$a==1 & d$t>=trt.time]) - mean(d$y[d$a==0 & d$t>=trt.time])

d = d[!(d$t %in% t_out),]

g6 = ggplot(data = d, aes(x = x, y = factor(t), fill = factor(a), color=factor(a)),fill="white") +
  stat_density_ridges(alpha=0.3, scale = .9,lwd=1.2, quantile_lines = TRUE, quantiles = 0.5) + 
  scale_fill_manual(values=c("0"=col1, "1"=col2), labels = c("C","T")) +
  scale_color_manual(values=c("0"=col1, "1"=col2), labels = c("C","T")) +
  theme_bw()+
  theme_ipsum_rc(plot_title_face = "bold") + #plain 
  xlab("X") + ylab("") + ggtitle("Scenario 6")+
  
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(axis.title.x = element_text(size=18),#margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=14),#margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=10),legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.background = element_rect(fill="white",colour ="white"),
        title = element_text(size=14))

library(cowplot)

pgrd <- plot_grid(
  g1 + theme(legend.position="none"),
  g2 + theme(legend.position="none"),
  g3 + theme(legend.position="none"),
  #align = 'vh',
  #hjust = -1,
  nrow = 1
)


plot_grid(pgrd, ncol = 1, rel_heights = c(1, .1))


pgrd <- plot_grid(
  g4 + theme(legend.position="none"),
  g5 + theme(legend.position="none"),
  g6 + theme(legend.position="none"),
  #align = 'vh',
  #hjust = -1,
  nrow = 1
)


plot_grid(pgrd, ncol = 1, rel_heights = c(1, .1))
