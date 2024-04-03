  
   # Differences in communication and incentives between high and low performers 
  
   - 2017 survey for students in test-taking grades. 
  
   ```{r, echo=FALSE, warning=FALSE, message=FALSE, fig.cap="Results for 4th Grade"} 
   library(stringr) 
  
   tab4th <- read.csv("https://raw.githubusercontent.com/maibennett/presentations/main/content/presentations/Attendance/Salem_20220405/images/survey_4th.csv") 
    
     tab10th <- read.csv("https://raw.githubusercontent.com/maibennett/presentations/main/content/presentations/Attendance/Salem_20220405/images/survey_10th.csv") 
      
       tab4th_long = tab4th %>% filter(GPA.Decile!="") %>% mutate(across(Told:Grades, ~as.numeric(str_replace_all(.x,"\\*","")))) 
        
         tab4th_long = tab4th_long %>% pivot_longer(!GPA.Decile, names_to = "label", values_to = "value") 
          
           tab4th_long = tab4th_long %>% group_by(label) %>% arrange(GPA.Decile, .by_group = TRUE) %>% 
               mutate(baseline = ifelse(GPA.Decile=="Baseline", value, NA), 
                                      baseline = mean(baseline, na.rm = TRUE), 
                                      value = ifelse(GPA.Decile == "Baseline", value, value + baseline)) %>% ungroup() 
            
             tab4th_long = tab4th_long %>% mutate(ci_upper = value + 0.004*1.67, 
                                                                                              ci_lower = value - 0.004*1.67) #most conservatives CIs given the SE 
            
             tab4th_long$label = factor(tab4th_long$label, levels = c("Grades", "Preparation", "Notification","Told")) 
              
               arrow_data = tab4th_long %>% group_by(label) %>%  
                   mutate(x = ifelse(GPA.Decile == "Baseline", value, NA), 
                                          x = mean(x, na.rm = TRUE), 
                                          xend = ifelse(GPA.Decile == "D1", value, NA), 
                                          xend = mean(xend, na.rm = TRUE), 
                                          xend = ifelse(x>xend, xend + 0.01, xend - 0.01)) %>% summarize(x = mean(x), 
                                                                                                                                                                        xend = mean(xend)) %>% ungroup() 
                
                
                 g1 = ggplot(data = tab4th_long %>% filter(GPA.Decile!="D10"), aes(y = label, x = value*100, color = GPA.Decile, 
                                                                                                                                                       fill = GPA.Decile, pch = GPA.Decile, size = GPA.Decile)) + 
                     geom_point() + xlab("% Response") + ylab("") + xlim(0,100)+ 
                     scale_color_manual(values = c("Baseline" = col_text, "D1" = cols[4], "D10" = cols[3]), name = "") + 
                     scale_fill_manual(values = c("Baseline" = col_text, "D1" = cols[4], "D10" = cols[3]), name = "") + 
                     scale_shape_manual(values = c("Baseline" = 22, "D1" = 21, "D10" = 21), name = "") + 
                     scale_size_manual(values = c("Baseline" = 2, "D1" = 4, "D10" = 4), name = "") + 
                         theme_ipsum_fsc() + 
                       theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"), 
                                              panel.grid.major.x = element_line(color = alpha("lightgrey",0.3)), 
                                              panel.grid.minor.x = element_blank(), 
                                              panel.grid.major.y = element_blank(), 
                                              panel.grid.minor.y = element_blank(), 
                                              axis.line = element_line(colour = "darkgrey") 
                                              )+ 
                       theme(legend.position = "right") + 
                       theme(axis.title.x = element_text(size=18),#margin = margin(t = 10, r = 0, b = 0, l = 0)), 
                                            axis.text.x = element_text(size=14), 
                                              axis.title.y = element_text(size=18),#margin = margin(t = 0, r = 10, b = 0, l = 0)), 
                                            axis.text.y = element_text(size=14), 
                                              plot.title = element_text(size=30), 
                                              legend.title = element_blank()) + 
                     geom_segment(data = arrow_data, aes(x = x*100, y = label, xend = xend*100, yend = label), 
                                                                arrow = arrow(type = "open", length = unit(0.1, "inches")), lineend = 'round', 
                                                        inherit.aes = FALSE, color = "darkgrey") 
                  
                  
                   x1 <- 80  # Example start point 
                   x2 <- 92  # Example end point 
                   y_level <- "Told"  # The y-level to place the rectangle 
                  
                   g1 + annotate("rect", xmin = x1, xmax = x2, ymin = 4 - 0.3, ymax = 4 + 0.3, alpha = 0, color = cols[2], linewidth = 1) 
                    
                     ``` 
                    
                    
                      
                     # Differences in communication and incentives between high and low performers 
                    
                     - 2017 survey for students in test-taking grades. 
                    
                     ```{r, echo=FALSE, warning=FALSE, message=FALSE, fig.cap="Results for 4th Grade"} 
                     arrow_data2 = tab4th_long %>% group_by(label) %>%  
                         mutate(x = ifelse(GPA.Decile == "Baseline", value, NA), 
                                                x = mean(x, na.rm = TRUE), 
                                                xend = ifelse(GPA.Decile == "D10", value, NA), 
                                                xend = mean(xend, na.rm = TRUE), 
                                                xend = ifelse(x>xend, xend + 0.01, xend - 0.01)) %>% summarize(x = mean(x), 
                                                                                                                                                                              xend = mean(xend)) %>% ungroup() 
                      
                      
                       g2 = ggplot(data = tab4th_long, aes(y = label, x = value*100, color = GPA.Decile, 
                                                                                                                               fill = GPA.Decile, pch = GPA.Decile, 
                                                                                                      size = GPA.Decile)) + 
                           geom_point() + xlab("% Response") + ylab("") + xlim(0,100)+ 
                           scale_color_manual(values = c("Baseline" = col_text, "D1" = cols[4], "D10" = cols[3]), name = "") + 
                           scale_fill_manual(values = c("Baseline" = col_text, "D1" = cols[4], "D10" = cols[3]), name = "") + 
                           scale_shape_manual(values = c("Baseline" = 22, "D1" = 21, "D10" = 21), name = "") + 
                           scale_size_manual(values = c("Baseline" = 2, "D1" = 4, "D10" = 4), name = "") + 
                               theme_ipsum_fsc() + 
                             theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"), 
                                                    panel.grid.major.x = element_line(color = alpha("lightgrey",0.3)), 
                                                    panel.grid.minor.x = element_blank(), 
                                                    panel.grid.major.y = element_blank(), 
                                                    panel.grid.minor.y = element_blank(), 
                                                    axis.line = element_line(colour = "darkgrey") 
                                                    )+ 
                             theme(legend.position = "right") + 
                             theme(axis.title.x = element_text(size=18),#margin = margin(t = 10, r = 0, b = 0, l = 0)), 
                                                  axis.text.x = element_text(size=14), 
                                                    axis.title.y = element_text(size=18),#margin = margin(t = 0, r = 10, b = 0, l = 0)), 
                                                  axis.text.y = element_text(size=14), 
                                                    plot.title = element_text(size=30), 
                                                    legend.title = element_blank()) + 
                           geom_segment(data = arrow_data, aes(x = x*100, y = label, xend = xend*100, yend = label), 
                                                                      arrow = arrow(type = "open", length = unit(0.1, "inches")), lineend = 'round', 
                                                              inherit.aes = FALSE, color = "darkgrey") + 
                             geom_segment(data = arrow_data2, aes(x = x*100, y = label, xend = xend*100, yend = label), 
                                                                        arrow = arrow(type = "open", length = unit(0.1, "inches")), lineend = 'round', 
                                                                inherit.aes = FALSE, color = "darkgrey") 
                        
                        
                         g2 
                        
                         ``` 