


carol <- read_excel(here("Data","Amostras_MICs_Completo.xlsx"))
carol$Ano2<-carol$Ano-2000
carol$Material<-factor(carol$Material,labels=c("Carcass","Food"))
carol$Projeto<-as.factor(carol$Projeto)


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


summary(mod1<-lmer(Ertapenem~Ano2:Sorovar:Material+Sorovar+Material+(Ano2|Projeto),data=carol))

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




#Predicted values
x5<-rep(0:15) #year


TC1<-mod1@beta[1]+mod1@beta[2]+(mod1@beta[5]*x5)
TF1<-mod1@beta[1]+mod1@beta[2]+mod1@beta[3]+(mod1@beta[7]*x5)
DC1<-mod1@beta[1]+(mod1@beta[4]*0:10)
DF1<-mod1@beta[1]+ mod1@beta[3]+(mod1@beta[6]*0:10)
DC1<-c(DC1,rep(NA,5))
DF1<-c(DF1,rep(NA,5))


pr<-cbind.data.frame(Year=x5,TC1,TF1,DC1,DF1)

pr1<-pr%>%
  gather(key="S:M",value="Halo",-Year)

pr1$Serovar<-ifelse(pr1$`S:M`=="TC1"|pr1$`S:M`=="TF1","Typhimurium", "Derby")
pr1$Material<-ifelse(pr1$`S:M`=="TC1"|pr1$`S:M`=="DC1","Carcass", "Food")


ggplot(pr1,aes(x=Year,y=Halo,col=Serovar))+
  theme_minimal()+
  theme(text=element_text(size=21))+
#  facet_grid(~Material)+
  geom_point(aes(shape = Material),size=3.5)+
  geom_line(aes(shape = Material),size=1.5)+
  scale_x_continuous(breaks = seq(0,16,2))+
  scale_y_continuous(breaks = seq(34,40,1))+
  xlab("Year")


#Multiple comparissons
betas1<-mod1@beta[c(-1,-2,-3)]
stderr1<-unlist(sqrt(diag(mod1@vcov_beta)))

ztest1<-cbind(betas1,stderr1[c(-1,-2,-3)])

z1<-2*pnorm(abs(ztest1[1]-ztest1[2])/sqrt(ztest1[1,2]^2+ztest1[2,2]^2),lower.tail = F)
2*pt(abs(ztest1[1]-ztest1[2])/sqrt(ztest1[1,2]^2+ztest1[2,2]^2),df=180+82-4,lower.tail = F)
  
z2<-2*pnorm(abs(ztest1[3]-ztest1[4])/sqrt(ztest1[3,2]^2+ztest1[4,2]^2),lower.tail = F)
2*pt(abs(ztest1[3]-ztest1[4])/sqrt(ztest1[3,2]^2+ztest1[4,2]^2),df=53+98-4,lower.tail = F)

z3<-2*pnorm(abs(ztest1[2]-ztest1[4])/sqrt(ztest1[2,2]^2+ztest1[4,2]^2),lower.tail = F)
2*pt(abs(ztest1[2]-ztest1[4])/sqrt(ztest1[2,2]^2+ztest1[4,2]^2),df=180+98-4,lower.tail = F)

z4<-2*pnorm(abs(ztest1[1]-ztest1[3])/sqrt(ztest1[1,2]^2+ztest1[3,2]^2),lower.tail = F)
2*pt(abs(ztest1[1]-ztest1[3])/sqrt(ztest1[1,2]^2+ztest1[3,2]^2),df=82+53-4,lower.tail = F)


summary(mod2<-lmer(Meropenem~Ano2:sorovar1:Material+sorovar1+Material+(Ano2|Projeto),data=carol))

Anova(mod2)

plot(fitted(mod2), resid(mod2))# this will create the plot
abline(0,0, col="red")

qqPlot(residuals(mod2))

shapiro.test(residuals(mod2))

TC2<-mod2@beta[1]+mod2@beta[2]+(mod2@beta[5]*x5)
TF2<-mod2@beta[1]+mod2@beta[2]+mod2@beta[3]+(mod2@beta[7]*x5)
DC2<-mod2@beta[1]+(mod2@beta[4]*0:10)
DF2<-mod2@beta[1]+ mod2@beta[3]+(mod2@beta[6]*0:10)
DC2<-c(DC2,rep(NA,5))
DF2<-c(DF2,rep(NA,5))




pr.1<-cbind.data.frame(Year=x5,TC2,TF2,DC2,DF2)

pr1.1<-pr.1%>%
  gather(key="S:M",value="Halo",-Year)

pr1.1$Serovar<-ifelse(pr1.1$`S:M`=="TC2"|pr1.1$`S:M`=="TF2","Typhimurium", "Derby")
pr1.1$Material<-ifelse(pr1.1$`S:M`=="TC2"|pr1.1$`S:M`=="DC2","Carcass", "Food")


ggplot(pr1.1,aes(x=Year,y=Halo,col=Serovar))+
  theme_minimal()+
  theme(text=element_text(size=21))+
  #  facet_grid(~Material)+
  geom_point(aes(shape = Material),size=3.5)+
  geom_line(aes(shape = Material),size=1.5)+
  scale_x_continuous(breaks = seq(0,16,2))+
  scale_y_continuous(breaks = seq(34,40,1))+
  xlab("Year")


#Multiple comparissons
betas2<-mod2@beta[c(-1,-2,-3)]
stderr2<-unlist(sqrt(diag(mod2@vcov_beta)))

ztest2<-cbind(betas2,stderr2[c(-1,-2,-3)])

z1.2<-2*pnorm(abs(ztest2[1]-ztest2[2])/sqrt(ztest2[1,2]^2+ztest2[2,2]^2),lower.tail = F)

z2.2<-2*pnorm(abs(ztest2[3]-ztest2[4])/sqrt(ztest2[3,2]^2+ztest2[4,2]^2),lower.tail = F)

z3.2<-2*pnorm(abs(ztest2[2]-ztest2[4])/sqrt(ztest2[2,2]^2+ztest2[4,2]^2),lower.tail = F)

z4.2<-2*pnorm(abs(ztest2[1]-ztest2[3])/sqrt(ztest2[1,2]^2+ztest2[3,2]^2),lower.tail = F)



summary(mod3<-lmer(Imipenem~Ano2:sorovar1:Material+sorovar1+Material+(Ano2|Projeto),data=carol))

Anova(mod3)

plot(fitted(mod3), resid(mod3))# this will create the plot
abline(0,0, col="red")

qqPlot(residuals(mod3))

shapiro.test(residuals(mod3))

TC3<-mod3@beta[1]+mod3@beta[2]+(mod3@beta[5]*x5)
TF3<-mod3@beta[1]+mod3@beta[2]+mod3@beta[3]+(mod3@beta[7]*x5)
DC3<-mod3@beta[1]+(mod3@beta[4]*0:10)
DF3<-mod3@beta[1]+ mod3@beta[3]+(mod3@beta[6]*0:10)
DC3<-c(DC3,rep(NA,5))
DF3<-c(DF3,rep(NA,5))


pr.2<-cbind.data.frame(Year=x5,TC3,TF3,DC3,DF3)

pr1.2<-pr.2%>%
  gather(key="S:M",value="Halo",-Year)

pr1.2$Serovar<-ifelse(pr1.2$`S:M`=="TC3"|pr1.2$`S:M`=="TF3","Typhimurium", "Derby")
pr1.2$Material<-ifelse(pr1.2$`S:M`=="TC3"|pr1.2$`S:M`=="DC3","Carcass", "Food")


ggplot(pr1.2,aes(x=Year,y=Halo,col=Serovar))+
  theme_minimal()+
  theme(text=element_text(size=21))+
  #  facet_grid(~Material)+
  geom_point(aes(shape = Material),size=3.5)+
  geom_line(aes(shape = Material),size=1.5)+
  scale_x_continuous(breaks = seq(0,16,2))+
  scale_y_continuous(breaks = seq(30,40,1))+
  xlab("Year")


#Multiple comparissons
betas3<-mod3@beta[c(-1,-2,-3)]
stderr3<-unlist(sqrt(diag(mod3@vcov_beta)))

ztest3<-cbind(betas3,stderr3[c(-1,-2,-3)])

z1.3<-2*pnorm(abs(ztest3[1]-ztest3[2])/sqrt(ztest3[1,2]^2+ztest3[2,2]^2),lower.tail = F) #Derby vs typh (carcass)

z2.3<-2*pnorm(abs(ztest3[3]-ztest3[4])/sqrt(ztest3[3,2]^2+ztest3[4,2]^2),lower.tail = F) #Derby vs typhi (food)

z3.3<-2*pnorm(abs(ztest3[2]-ztest3[4])/sqrt(ztest3[2,2]^2+ztest3[4,2]^2),lower.tail = F) # Carcass vs food (typhimu)

z4.3<-2*pnorm(abs(ztest3[1]-ztest3[3])/sqrt(ztest3[1,2]^2+ztest3[3,2]^2),lower.tail = F) # Carcass vs food (derby)





# All predictive plots together


colnames(pr)<-c("Year","Tc","TF","DC","DF")
colnames(pr.1)<-c("Year","Tc","TF","DC","DF")
colnames(pr.2)<-c("Year","Tc","TF","DC","DF")

pred_total<-rbind.data.frame(pr,pr.1,pr.2)
pred_total$ATM<-c(rep("Ertapenem",16),rep("Meropenem",16),rep("Imipenem",16))

pred_final<-pred_total%>%
  gather(key="Bacteria",value="Halo", -Year,-ATM)

pred_final$serovar<-ifelse(pred_final$Bacteria=="Tc"|pred_final$Bacteria=="TF","Typhimurium","Derby")
pred_final$matrix<-ifelse(pred_final$Bacteria=="Tc"|pred_final$Bacteria=="DC","Carcass","Food")

ggplot(pred_final,aes(x=Year,y=Halo,col=serovar))+
  theme_minimal()+
  theme(text=element_text(size=21))+
  facet_grid(matrix~ATM)+
  #geom_jitter(aes(shape = Material))+
  geom_jitter()+
  scale_x_continuous(breaks = seq(0,15,2))+
  xlab("Year")



