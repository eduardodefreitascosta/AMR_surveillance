


carol <- read_excel(here("Data","Amostras_MICs_Completo.xlsx"))
carol$Ano2<-carol$Ano-2000
carol$Material<-factor(carol$Material,labels=c("Carcass","Food"))
carol$Projeto<-as.factor(carol$Projeto)

carol<-subset(carol,ID != 296)

#################
# Carbapenem #
#################

# Histograms
ggplot(carol,aes(x=Ertapenem,fill=Sorovar))+
  theme_minimal()+
    geom_histogram("binwidth"=1,col="black",alpha=0.5)+
    facet_grid(~Material)

ggplot(carol,aes(x=Meropenem,fill=Sorovar))+
  theme_minimal()+
  geom_histogram("binwidth"=1,col="black",alpha=0.5)+
  facet_grid(~Material)

ggplot(carol,aes(x=Imipenem,fill=Sorovar))+
  theme_minimal()+
  geom_histogram("binwidth"=1,col="black",alpha=0.5)+
  facet_grid(~Material)


#Scatter-plots
ggplot(carol,aes(x=Ano,y=Ertapenem,col=Sorovar))+
  theme_minimal()+
  facet_grid(~Material)+
  geom_jitter()+
  geom_smooth(method = "lm")+
  scale_x_continuous(breaks = seq(2000,2015,2))

ggplot(carol,aes(x=Ano,y=Meropenem,col=Sorovar))+
  theme_minimal()+
  facet_grid(~Material)+
  geom_jitter()+
  geom_smooth(method = "lm")+
  scale_x_continuous(breaks = seq(2000,2015,2))

ggplot(carol,aes(x=Ano,y=Imipenem,col=Sorovar))+
  theme_minimal()+
  facet_grid(~Material)+
  #geom_jitter(aes(shape = Material))+
  geom_jitter()+
  geom_smooth(method = "lm")+
  scale_x_continuous(breaks = seq(2000,2015,2))


geral<-carol[,c(3:11,16)]%>%
  gather(key="Bacteria",value="Halo", -Ano2,-Sorovar,-Projeto,-`Trabalho de origem`,-Responsável,-Paper,-Material)


ggplot(geral,aes(x=Ano2,y=Halo,col=Sorovar))+
  theme_minimal()+
  theme(text=element_text(size=21))+
  facet_grid(Material~Bacteria)+
  #geom_jitter(aes(shape = Material))+
  geom_jitter()+
  geom_smooth(method = "lm",size=2)+
  scale_x_continuous(breaks = seq(0,15,2))+
  xlab("Year")


#################
# Linear models #
#################


### Ertapenem ###


summary(mod1<-lmer(Ertapenem~Ano2:Sorovar:Material+Sorovar+Material+(1|Projeto),data=carol))

#ICC(mod1) Using SAS PROC MIXED to fit multilevel models, hierarchical models, and individual growth models
# https://www.biostat.jhsph.edu/~fdominic/teaching/bio656/lectures/4.linear.randomcoeff.pdf
# http://www.bristol.ac.uk/cmm/learning/videos/random-slopes.html


anova(mod1)

plot(fitted(mod1), resid(mod1))# this will create the plot
abline(0,0, col="red")

qqPlot(residuals(mod1))
acf(residuals(mod1))
hist(residuals(mod1))
shapiro.test(residuals(mod1))





### Meropenem ###

summary(mod2<-lmer(Meropenem~Ano2:Sorovar:Material+Sorovar+Material+(1|Projeto),data=carol))

Anova(mod2)

plot(fitted(mod2), resid(mod2))# this will create the plot
abline(0,0, col="red")

qqPlot(residuals(mod2))

shapiro.test(residuals(mod2))


### Imipenem ###
summary(mod3<-lmer(Imipenem~Ano2:Sorovar:Material+Sorovar+Material+(1|Projeto),data=carol))

Anova(mod3)

plot(fitted(mod3), resid(mod3))# this will create the plot
abline(0,0, col="red")

qqPlot(residuals(mod3))

shapiro.test(residuals(mod3))



# All predictive plots together


carol$pred.erta<-predict(mod1, re.form = NA)
carol$pred.mero<-predict(mod2, re.form = NA)
carol$pred.imi<-predict(mod3, re.form = NA)

pred_final1<-carol[,c(3,8,9:11,16)]%>%
  gather(key="ATM",value="Halo",-Sorovar,-Material,-Ano2)

pred_final2<-carol[,c(3,8,17:19,16)]%>%
  gather(key="ATM",value="pred",-Sorovar,-Material,-Ano2)

pred_final2$ATM<-ifelse(pred_final2$ATM=="pred.erta","Ertapenem",
                        ifelse(pred_final2$ATM=="pred.imi","Imipenem","Meropenem"))

pred_final<-left_join(pred_final1,pred_final2)



mat.labs <- c("Animal", "Food")
names(mat.labs) <- c("Carcass", "Food")


ggplot(pred_final1,aes(x=Ano2,y=Halo,col=Sorovar))+
  theme_minimal()+
  theme(text=element_text(size=20))+
  facet_grid(Material~ATM,scales = "free_y",labeller=labeller(Material=mat.labs) )+
  geom_jitter()+
  geom_line(data = pred_final2, aes(x = Ano2, y = pred),size=1)+
  scale_x_continuous(breaks = seq(0,15,2))+
  labs(x ="Year", y = "Halo diameter (mm)")+
  scale_color_discrete(name = "Serovar", labels = c("Derby", "Typhimurium"))+
  theme(panel.spacing.x = unit(1, "lines"))+
  theme(panel.spacing.y = unit(1.5, "lines"))
ggsave(here("Figures",'Halo.png'), dpi = 300, height = 10, width = 15, unit = 'in',bg="white")




ggplot(pred_final1,aes(x=Ano2,y=Halo,col=Sorovar))+
  theme_minimal()+
  theme(text=element_text(size=20))+
  facet_grid(ATM~Material,scales = "free_y",labeller=labeller(Material=mat.labs) )+
  geom_jitter()+
  geom_line(data = pred_final2, aes(x = Ano2, y = pred),size=1)+
  scale_x_continuous(breaks = seq(0,15,2))+
  labs(x ="Year", y = "Halo diameter (mm)")+
  scale_color_discrete(name = "Serovar", labels = c("Derby", "Typhimurium"))+
  theme(panel.spacing.x = unit(1, "lines"))+
  theme(panel.spacing.y = unit(1.5, "lines"))
ggsave(here("Figures",'Halo_1.png'), dpi = 300, height = 10, width = 15, unit = 'in',bg="white")





#Export tables

s.mod1<-summary(mod1)
s.mod2<-summary(mod2)
s.mod3<-summary(mod3)
  

final_res<-list(
  erta=round(cbind(s.mod1$coefficients[,c(1,5)],confint(mod1)[3:9,]),2),
  erta_var=s.mod1$varcor,
  mero=round(cbind(s.mod2$coefficients[,c(1,5)],confint(mod2)[3:9,]),2),
  mer_var=s.mod2$varcor,
  imi=round(cbind(s.mod3$coefficients[,c(1,5)],confint(mod3)[3:9,]),2),
  imi_var=s.mod3$varcor
)


sink(here("Outputs","Halo_reg.txt"))
print(final_res)
sink() 
