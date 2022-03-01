

# Creating the dataset
carol<-read_excel(here("Data","MIC_reg.xlsx"),sheet = "Todas")

carol$sorovar<-factor(carol$sorovar,labels=c("Typh","Derby"))
carol$material<-factor(carol$material,labels=c("Carc","Food"))

#Predicted values
x5<-rep(0:15) #year


###############
# Descriptive #
###############


medias<-carol%>%
   group_by(ano)%>%
   summarise(med_coli=mean(med_coli),med_cipro=mean(med_cipro),med_cefa=mean(med_ceftaz),med_cefat=mean(med_cefo))


medias1<-medias%>%
   gather(key="ATB",value="MIC_mean",-ano)


m1<-ggplot(medias1,aes(x=ano,y=MIC_mean,color=ATB))+
   theme_minimal()+
   geom_jitter()+
   geom_smooth(method = "lm")
   


MIC<-carol[,c(3,4,6,14,18,22,26)]

MIC1<-MIC%>%
   gather(key="ATB",value="MIC",-ano,-sorovar,-material)


m1<-ggplot(carol,aes(x=int_coli))+
   theme_minimal()+
   geom_bar()+
   xlab("MIC interval")+
   ylab("")+
   ggtitle("Colistina")

m2<-ggplot(carol,aes(x=int_cipro))+
   theme_minimal()+
   geom_bar()+
   xlab("MIC interval")+
   ylab("")+
   ggtitle("Ciprofloxacin")

m3<-ggplot(carol,aes(x=int_ceftaz))+
   theme_minimal()+
   geom_bar()+
   xlab("MIC interval")+
   ylab("")+
   ggtitle("Cefatazidime")

m4<-ggplot(carol,aes(x=int_cefo))+
   theme_minimal()+
   geom_bar()+
   xlab("MIC interval")+
   ylab("")+
   ggtitle("Cefotaxime")

grid.arrange(m1,m2,m3,m4)




#####################
#Determinsitic model#
#####################

# A Study of Interval Censoring in Parametric Regression Models. 10.1023/A:1009681919084


mod1<-survreg(Surv(time=min_coli,time2=max_coli,type="interval2")~ 
                 ano:sorovar:material+sorovar+material,data=carol,cluster=projeto,dist="exponential")
summary(mod1)

linearHypothesis(mod1, "ano") #TC
linearHypothesis(mod1, "ano+ano:materialFood") #TF
linearHypothesis(mod1, "ano+ano:sorovarDerby") #DC
linearHypothesis(mod1, "ano+ano:materialFood+ano:sorovarDerby+ano:sorovarDerby:materialFood") #DF

TC1<-exp(mod1$coefficients[1]+mod1$coefficients[2]*x5) 

TF1<-exp(mod1$coefficients[1]+mod1$coefficients[2]*x5+mod1$coefficients[4]+mod1$coefficients[6]*x5)

DC1<-exp(mod1$coefficients[1]+mod1$coefficients[3]+(mod1$coefficients[2]+mod1$coefficients[5])*0:10)

DF1<-exp(mod1$coefficients[1]+mod1$coefficients[3]+mod1$coefficients[4]+mod1$coefficients[7]+
            (mod1$coefficients[2]+mod1$coefficients[5]+mod1$coefficients[6]+mod1$coefficients[8])*0:10)


DC1<-c(DC1,rep(NA,5))
DF1<-c(DF1,rep(NA,5))

pr.1<-cbind.data.frame(Year=x5,TC=TC1,TF=TF1,DC=DC1,DF=DF1)

pr1<-pr.1%>%
   gather(key="S:M",value="Mean MIC",-Year)

pr1$Serovar<-ifelse(pr1$`S:M`=="TC1"|pr1$`S:M`=="TF1","Typhimurium", "Derby")
pr1$Material<-ifelse(pr1$`S:M`=="TC1"|pr1$`S:M`=="DC1","Carcass", "Food")


ggplot(pr1,aes(x=Year,y=`Mean MIC`,col=Serovar))+
   theme_minimal()+
   theme(text=element_text(size=21))+
   #  facet_grid(~Material)+
   geom_point(aes(shape = Material),size=3.5)+
   geom_line(aes(shape = Material),size=1.5)+
   scale_x_continuous(breaks = seq(0,16,2))+
   scale_y_continuous(breaks = seq(0,4,0.5))+
   xlab("Year")

ggsave(here("Figures","Picture6.png"))


carol$max_cipro1<-ifelse(carol$max_cipro=="Inf"," ",as.numeric(carol$max_cipro))
carol$max_cipro1<-as.numeric(carol$max_cipro1)

mod2<-survreg(Surv(time=min_cipro,time2=max_cipro1,type="interval2")~
                 ano:sorovar:material+sorovar+material,data=carol,cluster=projeto,dist="exponential")
summary(mod2)


linearHypothesis(mod2, "ano") #TC
linearHypothesis(mod2, "ano+ano:materialFood") #TF
linearHypothesis(mod2, "ano+ano:sorovarDerby") #DC
linearHypothesis(mod2, "ano+ano:materialFood+ano:sorovarDerby+ano:sorovarDerby:materialFood") #DF


TC2<-exp(mod2$coefficients[1]+mod2$coefficients[2]*x5) 

TF2<-exp(mod2$coefficients[1]+mod2$coefficients[4]+(mod2$coefficients[2]+mod2$coefficients[6])*x5)

DC2<-exp(mod2$coefficients[1]+mod2$coefficients[3]+(mod2$coefficients[2]+mod2$coefficients[5])*0:10)

DF2<-exp(mod2$coefficients[1]+mod2$coefficients[3]+mod2$coefficients[4]+mod2$coefficients[7]+
            (mod2$coefficients[2]+mod2$coefficients[5]+mod2$coefficients[6]+mod2$coefficients[8])*0:10)

DC2<-c(DC2,rep(NA,5))
DF2<-c(DF2,rep(NA,5))



pr.2<-cbind.data.frame(Year=x5,TC=TC2,TF=TF2,DC=DC2,DF=DF2)

pr2<-pr.2%>%
   gather(key="S:M",value="Mean MIC",-Year)

pr2$Serovar<-ifelse(pr2$`S:M`=="TC2"|pr2$`S:M`=="TF2","Typhimurium", "Derby")
pr2$Material<-ifelse(pr2$`S:M`=="TC2"|pr2$`S:M`=="DC2","Carcass", "Food")


ggplot(pr2,aes(x=Year,y=`Mean MIC`,col=Serovar))+
   theme_minimal()+
   theme(text=element_text(size=21))+
   #  facet_grid(~Material)+
   geom_point(aes(shape = Material),size=3.5)+
   geom_line(aes(shape = Material),size=1.5)+
   scale_x_continuous(breaks = seq(0,16,2))+
   scale_y_continuous(breaks = seq(0,1,0.05))+
   xlab("Year")

ggsave(here("Figures","Picture7.png"))




mod3<-survreg(Surv(time=min_ceftaz,time2=max_ceftaz,type="interval2")~ano:sorovar:material+sorovar+material,cluster=projeto,data=carol,dist="exponential")
summary(mod3)

linearHypothesis(mod3, "ano") #TC
linearHypothesis(mod3, "ano+ano:materialFood") #TF
linearHypothesis(mod3, "ano+ano:sorovarDerby") #DC
linearHypothesis(mod3, "ano+ano:materialFood+ano:sorovarDerby+ano:sorovarDerby:materialFood") #DF

TC3<-exp(mod3$coefficients[1]+mod3$coefficients[2]*x5) 

TF3<-exp(mod3$coefficients[1]+mod3$coefficients[4]+(mod3$coefficients[2]+mod3$coefficients[6])*x5)

DC3<-exp(mod3$coefficients[1]+mod3$coefficients[3]+(mod3$coefficients[2]+mod3$coefficients[5])*0:10)

DF3<-exp(mod3$coefficients[1]+mod3$coefficients[3]+mod3$coefficients[4]+mod3$coefficients[7]+
            (mod3$coefficients[2]+mod3$coefficients[5]+mod3$coefficients[6]+mod3$coefficients[8])*0:10)
DC3<-c(DC3,rep(NA,5))
DF3<-c(DF3,rep(NA,5))


pr.3<-cbind.data.frame(Year=x5,TC=TC3,TF=TF3,DC=DC3,DF=DF3)

pr3<-pr.3%>%
   gather(key="S:M",value="Mean MIC",-Year)

pr3$Serovar<-ifelse(pr3$`S:M`=="TC3"|pr3$`S:M`=="TF3","Typhimurium", "Derby")
pr3$Material<-ifelse(pr3$`S:M`=="TC3"|pr3$`S:M`=="DC3","Carcass", "Food")


ggplot(pr3,aes(x=Year,y=`Mean MIC`,col=Serovar))+
   theme_minimal()+
   theme(text=element_text(size=21))+
   #  facet_grid(~Material)+
   geom_point(aes(shape = Material),size=3.5)+
   geom_line(aes(shape = Material),size=1.5)+
   scale_x_continuous(breaks = seq(0,16,2))+
   scale_y_continuous(breaks = seq(0,1,0.05))+
   xlab("Year")

ggsave(here("Figures","Picture8.png"))





mod4<-survreg(Surv(time=min_cefo,time2=max_cefo,type="interval2")~ano:sorovar:material+sorovar+material,cluster=projeto,data=carol,dist="exponential")
summary(mod4)

linearHypothesis(mod4, "ano") #TC
linearHypothesis(mod4, "ano+ano:materialFood") #TF
linearHypothesis(mod4, "ano+ano:sorovarDerby") #DC
linearHypothesis(mod4, "ano+ano:materialFood+ano:sorovarDerby+ano:sorovarDerby:materialFood") #DF

TC4<-exp(mod4$coefficients[1]+mod4$coefficients[2]*x5) 

TF4<-exp(mod4$coefficients[1]+mod4$coefficients[4]+(mod4$coefficients[2]+mod4$coefficients[6])*x5)

DC4<-exp(mod4$coefficients[1]+mod4$coefficients[3]+(mod4$coefficients[2]+mod4$coefficients[5])*0:10)

DF4<-exp(mod4$coefficients[1]+mod4$coefficients[3]+mod4$coefficients[4]+mod4$coefficients[7]+
            (mod4$coefficients[2]+mod4$coefficients[5]+mod4$coefficients[6]+mod4$coefficients[8])*0:10)

DC4<-c(DC4,rep(NA,5))
DF4<-c(DF4,rep(NA,5))



pr.4<-cbind.data.frame(Year=x5,TC=TC4,TF=TF4,DC=DC4,DF=DF4)

pr4<-pr.4%>%
   gather(key="S:M",value="Mean MIC",-Year)

pr4$Serovar<-ifelse(pr4$`S:M`=="TC4"|pr4$`S:M`=="TF4","Typhimurium", "Derby")
pr4$Material<-ifelse(pr4$`S:M`=="TC4"|pr4$`S:M`=="DC4","Carcass", "Food")


ggplot(pr4,aes(x=Year,y=`Mean MIC`,col=Serovar))+
   theme_minimal()+
   theme(text=element_text(size=21))+
   #  facet_grid(~Material)+
   geom_point(aes(shape = Material),size=3.5)+
   geom_line(aes(shape = Material),size=1.5)+
   scale_x_continuous(breaks = seq(0,16,2))+
   scale_y_continuous(breaks = seq(0,0.2,0.025))+
   xlab("Year")
ggsave(here("Figures","Picture9.png"))


#All plots together

pred<-rbind.data.frame(pr.1,pr.2,pr.3,pr.4)
pred$ATM<-c(rep("Colistin",16),rep("Ciprofloxacin",16),rep("Cefatazidim",16),rep("Cefotaxime",16))

pred_total<-pred%>%
   gather(key="Bacteria",value="MIC", -Year,-ATM)

pred_total$serovar<-ifelse(pred_total$Bacteria=="TC"|pred_total$Bacteria=="TF","Typh","Derby")
pred_total$matrix<-ifelse(pred_total$Bacteria=="TC"|pred_total$Bacteria=="DC","Carc","Food")


carol1<-carol[,c(3,4,6,13,17,21,25)]

carol2<-carol1%>%
   gather(key="ATM1",value="Med_MIC",-ano,-sorovar,-material)

colnames(carol2)[1]<-"Year"
colnames(carol2)[3]<-"matrix"
colnames(carol2)[2]<-"serovar"


carol2$ATM<-ifelse(carol2$ATM1=="med_coli","Colistin",
                   ifelse(carol2$ATM1=="med_cipro","Ciprofloxacin",
                          ifelse(carol2$ATM1=="med_cefo","Cefotaxime","Cefatazidim")))


geral_final<-left_join(carol2,pred_total,by=c("Year","ATM","serovar","matrix"))


ggplot(geral_final,aes(x=Year,y=Med_MIC,col=serovar))+
   theme_minimal()+
   theme(text=element_text(size=21))+
   facet_grid(matrix~ATM)+
   geom_jitter()+
   geom_line(aes(x=Year,y = MIC), size = 1)+
   scale_x_continuous(breaks = seq(0,16,2))+
   xlab("Year")

