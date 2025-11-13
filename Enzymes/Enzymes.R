enz <- read.csv("~/Desktop/Enzymes.csv")
str(enz)


require(tidyverse)

ggplot(enz %>% filter(Light == 8 | Day == 0)) +
  geom_point(aes(Day, Glucosidase/Aminopeptidase, colour = factor(Temperature))) +
  geom_smooth(aes(Day, Glucosidase/Aminopeptidase, colour = factor(Temperature)))



require(mgcv)
# ?mgcv
# 
# m1 <- gamm(Glucosidase ~ Day * Temperature, random = ~1|Tank, data = enz[!enz$Light == "Dark",])
# 
# m1 <- gam(Glucosidase+0.1 ~ s(Day, Temperature, k = 5, bs = "fs", m = 2), family = Gamma(link = "log"),
#           data = enz[!enz$Light == "Dark",], method = "REML")
# hist(resid(m1))


mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .3, .2, .5),"cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15, hjust = 0),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.key.height = unit(.45, "cm"),
                 legend.key.width = unit(.45, "cm"),
                 legend.spacing.x = unit(.1, "cm"),
                 legend.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title.align = 0,
                 legend.title = element_blank(),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 15, hjust = 0),
                 panel.spacing = unit(.5, "cm"),
                 text = element_text(family = "Helvetica Neue"))

Gstats <- with(enz %>% filter(Light == 8 | Day == 0), 
              data.frame(Day = aggregate(Day, by = list(Group), mean)[,2],
                         Temperature = aggregate(Temperature, by = list(Group), mean)[,2],
                         Group = aggregate(Glucosidase, by = list(Group), mean)[,1],
                         Mean = aggregate(Glucosidase, by = list(Group), mean)[,2],
                         SD = aggregate(Glucosidase, by = list(Group), sd)[,2],
                         N = aggregate(Glucosidase, by = list(Group), length)[,2]))
Gstats$SEM <- with(Gstats, SD/sqrt(N))

Gstats2 <- with(enz %>% filter(Light == 8 | Day == 0), 
              data.frame(Day = aggregate(Glucosidase, by = list(Day, Temperature), mean)[,1],
                         Temperature = aggregate(Glucosidase, by = list(Day, Temperature), mean)[,2],
                         Mean = aggregate(Glucosidase, by = list(Day, Temperature), mean)[,3],
                         SD = aggregate(Glucosidase, by = list(Day, Temperature), sd)[,3],
                         N = aggregate(Glucosidase, by = list(Day, Temperature), length)[,3]))
Gstats2$Nstat <- with(Gstats2, N/3)
Gstats2$SEM <- with(Gstats2, SD/sqrt(Nstat))
Gstats2 <- rbind(with(enz %>% filter(Day == 0),
                     data.frame(Day = 0, Temperature = NA, Mean = mean(Glucosidase),
                                SD = sd(Glucosidase), N = length(Glucosidase), Nstat = length(Glucosidase)/3,
                                SEM = sd(Glucosidase)/sqrt(length(Glucosidase)/3))),
                Gstats2)

Gp <- ggplot() +
  geom_point(data = enz %>% filter(Light == 8 | Day == 0),
             aes(Day, Glucosidase, group = factor(Group), colour = factor(Temperature)),
             position = position_dodge(width = 6), shape = 16, size = 3, alpha = 0.2) +
  geom_pointrange(data = Gstats,
                  aes(Day, Mean, ymin = Mean - SEM, ymax = Mean + SEM,
                      group = factor(Group), colour = factor(Temperature)),
                  position = position_dodge(width = 6),
                  alpha = 0.5, shape = 16, size = 0.7) +
  geom_pointrange(data = Gstats2,
                  aes(Day, Mean, ymin = Mean - SEM, ymax = Mean + SEM,
                      colour = factor(Temperature)), size = 1) +
  geom_smooth(data = enz %>% filter(Light == 8 | Day == 0),
              aes(Day, Glucosidase, colour = factor(Temperature)), se = F, span = 0.5) +
  scale_colour_manual(breaks = c(15, 20, 25), na.value = "#000000",
                      values = c("#6ea4be", "#f5a54a", "#b3061e"),
                      labels = c("15°C", "20°C", "25°C")) +
  labs(x = "Detrital age (d)", 
       y = expression(beta*"-glucosidase activity ("*mu*"g cm"^-2*" h"^-1*")")) +
  coord_cartesian(ylim = c(0, 8), xlim = c(-5, 120), expand = F, clip = "off") +
  mytheme + theme(legend.position = c(0.9, 0.9))
  

Astats <- with(enz %>% filter(Light == 8 | Day == 0), 
              data.frame(Day = aggregate(Day, by = list(Group), mean)[,2],
                         Temperature = aggregate(Temperature, by = list(Group), mean)[,2],
                         Group = aggregate(Aminopeptidase, by = list(Group), mean)[,1],
                         Mean = aggregate(Aminopeptidase, by = list(Group), mean)[,2],
                         SD = aggregate(Aminopeptidase, by = list(Group), sd)[,2],
                         N = aggregate(Aminopeptidase, by = list(Group), length)[,2]))
Astats$SEM <- with(Astats, SD/sqrt(N))

Astats2 <- with(enz %>% filter(Light == 8 | Day == 0), 
               data.frame(Day = aggregate(Aminopeptidase, by = list(Day, Temperature), mean)[,1],
                          Temperature = aggregate(Aminopeptidase, by = list(Day, Temperature), mean)[,2],
                          Mean = aggregate(Aminopeptidase, by = list(Day, Temperature), mean)[,3],
                          SD = aggregate(Aminopeptidase, by = list(Day, Temperature), sd)[,3],
                          N = aggregate(Aminopeptidase, by = list(Day, Temperature), length)[,3]))
Astats2$Nstat <- with(Astats2, N/3)
Astats2$SEM <- with(Astats2, SD/sqrt(Nstat))
Astats2 <- rbind(with(enz %>% filter(Day == 0),
                     data.frame(Day = 0, Temperature = NA, Mean = mean(Aminopeptidase),
                                SD = sd(Aminopeptidase), N = length(Aminopeptidase), Nstat = length(Aminopeptidase)/3,
                                SEM = sd(Aminopeptidase)/sqrt(length(Aminopeptidase)/3))),
                Astats2)

Ap <- ggplot() +
  geom_point(data = enz %>% filter(Light == 8 | Day == 0),
             aes(Day, Aminopeptidase, group = factor(Group), colour = factor(Temperature)),
             position = position_dodge(width = 6), shape = 16, size = 3, alpha = 0.2) +
  geom_pointrange(data = Astats,
                  aes(Day, Mean, ymin = Mean - SEM, ymax = Mean + SEM,
                      group = factor(Group), colour = factor(Temperature)),
                  position = position_dodge(width = 6),
                  alpha = 0.5, shape = 16, size = 0.7) +
  geom_pointrange(data = Astats2,
                  aes(Day, Mean, ymin = Mean - SEM, ymax = Mean + SEM,
                      colour = factor(Temperature)), size = 1) +
  geom_smooth(data = enz %>% filter(Light == 8 | Day == 0),
              aes(Day, Aminopeptidase, colour = factor(Temperature)), se = F, span = 0.44) +
  scale_colour_manual(breaks = c(15, 20, 25), na.value = "#000000",
                      values = c("#6ea4be", "#f5a54a", "#b3061e"),
                      labels = c("15°C", "20°C", "25°C")) +
  labs(x = "Detrital age (d)", 
       y = expression("Aminopeptidase activity ("*mu*"g cm"^-2*" h"^-1*")")) +
  coord_cartesian(ylim = c(0, 100), xlim = c(-5, 120), expand = F, clip = "off") +
  mytheme + theme(legend.position = "none")

require(cowplot)
plot_grid(Gp + theme(axis.title.x = element_blank(), 
                     axis.text.x = element_blank()), 
          Ap, nrow = 2, labels = "auto", label_fontfamily = "Helvetica Neue", 
          label_size = 15, align = "v",
          rel_heights = c(0.931, 1))



LGstats <- with(enz %>% filter(Temperature == 15 | Day == 0), 
               data.frame(Day = aggregate(Day, by = list(Group), mean)[,2],
                          Light = aggregate(Light, by = list(Group), mean)[,2],
                          Group = aggregate(Glucosidase, by = list(Group), mean)[,1],
                          Mean = aggregate(Glucosidase, by = list(Group), mean)[,2],
                          SD = aggregate(Glucosidase, by = list(Group), sd)[,2],
                          N = aggregate(Glucosidase, by = list(Group), length)[,2]))
LGstats$SEM <- with(LGstats, SD/sqrt(N))

LGstats2 <- with(enz %>% filter(Temperature == 15 | Day == 0), 
                data.frame(Day = aggregate(Glucosidase, by = list(Day, Light), mean)[,1],
                           Light = aggregate(Glucosidase, by = list(Day, Light), mean)[,2],
                           Mean = aggregate(Glucosidase, by = list(Day, Light), mean)[,3],
                           SD = aggregate(Glucosidase, by = list(Day, Light), sd)[,3],
                           N = aggregate(Glucosidase, by = list(Day, Light), length)[,3]))
LGstats2$Nstat <- with(LGstats2, N/3)
LGstats2$SEM <- with(LGstats2, SD/sqrt(Nstat))
LGstats2 <- rbind(with(enz %>% filter(Day == 0),
                      data.frame(Day = 0, Light = NA, Mean = mean(Glucosidase),
                                 SD = sd(Glucosidase), N = length(Glucosidase), Nstat = length(Glucosidase)/3,
                                 SEM = sd(Glucosidase)/sqrt(length(Glucosidase)/3))),
                  LGstats2)

LGp <- ggplot() +
  geom_point(data = enz %>% filter(Temperature == 15 | Day == 0),
             aes(Day, Glucosidase, group = factor(Group), colour = factor(Light)),
             position = position_dodge(width = 6), shape = 16, size = 3, alpha = 0.2) +
  geom_pointrange(data = LGstats,
                  aes(Day, Mean, ymin = Mean - SEM, ymax = Mean + SEM,
                      group = factor(Group), colour = factor(Light)),
                  position = position_dodge(width = 6),
                  alpha = 0.5, shape = 16, size = 0.7) +
  geom_pointrange(data = LGstats2,
                  aes(Day, Mean, ymin = Mean - SEM, ymax = Mean + SEM,
                      colour = factor(Light)), size = 1) +
  geom_smooth(data = enz %>% filter(Temperature == 15 | Day == 0),
              aes(Day, Glucosidase, colour = factor(Light)), se = F, span = 0.4) +
  scale_colour_manual(breaks = c(0, 8), na.value = "#000000",
                      values = c("#6ea4be", "#f5a54a"),
                      labels = c("0 µmol photons m–2 s–1",
                                 "8 µmol photons m–2 s–1")) +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  labs(x = "Detrital age (d)", 
       y = expression(beta*"-glucosidase activity ("*mu*"g cm"^-2*" h"^-1*")")) +
  coord_cartesian(ylim = c(0, 10), xlim = c(-5, 120), expand = F, clip = "off") +
  mytheme + theme(legend.position = c(0.8, 0.9))


LAstats <- with(enz %>% filter(Temperature == 15 | Day == 0), 
                data.frame(Day = aggregate(Day, by = list(Group), mean)[,2],
                           Light = aggregate(Light, by = list(Group), mean)[,2],
                           Group = aggregate(Aminopeptidase, by = list(Group), mean)[,1],
                           Mean = aggregate(Aminopeptidase, by = list(Group), mean)[,2],
                           SD = aggregate(Aminopeptidase, by = list(Group), sd)[,2],
                           N = aggregate(Aminopeptidase, by = list(Group), length)[,2]))
LAstats$SEM <- with(LAstats, SD/sqrt(N))

LAstats2 <- with(enz %>% filter(Temperature == 15 | Day == 0), 
                 data.frame(Day = aggregate(Aminopeptidase, by = list(Day, Light), mean)[,1],
                            Light = aggregate(Aminopeptidase, by = list(Day, Light), mean)[,2],
                            Mean = aggregate(Aminopeptidase, by = list(Day, Light), mean)[,3],
                            SD = aggregate(Aminopeptidase, by = list(Day, Light), sd)[,3],
                            N = aggregate(Aminopeptidase, by = list(Day, Light), length)[,3]))
LAstats2$Nstat <- with(LAstats2, N/3)
LAstats2$SEM <- with(LAstats2, SD/sqrt(Nstat))
LAstats2 <- rbind(with(enz %>% filter(Day == 0),
                       data.frame(Day = 0, Light = NA, Mean = mean(Aminopeptidase),
                                  SD = sd(Aminopeptidase), N = length(Aminopeptidase), Nstat = length(Aminopeptidase)/3,
                                  SEM = sd(Aminopeptidase)/sqrt(length(Aminopeptidase)/3))),
                  LAstats2)

LAp <- ggplot() +
  geom_point(data = enz %>% filter(Temperature == 15 | Day == 0),
             aes(Day, Aminopeptidase, group = factor(Group), colour = factor(Light)),
             position = position_dodge(width = 6), shape = 16, size = 3, alpha = 0.2) +
  geom_pointrange(data = LAstats,
                  aes(Day, Mean, ymin = Mean - SEM, ymax = Mean + SEM,
                      group = factor(Group), colour = factor(Light)),
                  position = position_dodge(width = 6),
                  alpha = 0.5, shape = 16, size = 0.7) +
  geom_pointrange(data = LAstats2,
                  aes(Day, Mean, ymin = Mean - SEM, ymax = Mean + SEM,
                      colour = factor(Light)), size = 1) +
  geom_smooth(data = enz %>% filter(Temperature == 15 | Day == 0),
              aes(Day, Aminopeptidase, colour = factor(Light)), se = F, span = 0.3) +
  scale_colour_manual(breaks = c(0, 8), na.value = "#000000",
                      values = c("#6ea4be", "#f5a54a"),
                      labels = c("0 µmol photons m–2 s–1",
                                 "8 µmol photons m–2 s–1")) +
  labs(x = "Detrital age (d)", 
       y = expression("Aminopeptidase activity ("*mu*"g cm"^-2*" h"^-1*")")) +
  coord_cartesian(ylim = c(0, 50), xlim = c(-5, 120), expand = F, clip = "off") +
  mytheme + theme(legend.position = "none")

plot_grid(LGp + theme(axis.title.x = element_blank(), 
                     axis.text.x = element_blank()), 
          LAp, nrow = 2, labels = "auto", label_fontfamily = "Helvetica Neue", 
          label_size = 15, align = "v",
          rel_heights = c(0.931, 1))

