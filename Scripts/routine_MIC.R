

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



############
# Colistin #
############

carol_coli<-subset(carol,id != 296 & id != 281 )


mod1<-survreg(Surv(time=min_coli,time2=max_coli,type="interval2")~ 
                 ano:sorovar:material+sorovar+material,data=carol_coli,cluster = projeto,dist = "exponential")



s.mod1<-summary(mod1)

carol_coli$pred<-predict(mod1)

ggplot(carol_coli,aes(y=med_coli,x=ano1,color=sorovar))+
   facet_grid(~material)+
   geom_jitter()+
   geom_line(data = carol_coli, aes(x = ano1, y = pred))
   



#################
# Ciprofoxacina #
################


carol_cipro<-subset(carol,id != 296 & id != 8 & id != 254)
carol_cipro$max_cipro<-ifelse(carol_cipro$max_cipro=="Inf"," ",as.numeric(carol$max_cipro))
carol_cipro$max_cipro<-as.numeric(carol_cipro$max_cipro)

mod2<-survreg(Surv(time=min_cipro,time2=max_cipro,type="interval2")~
                 ano:sorovar:material+sorovar+material,data=carol_cipro,
              cluster=projeto,dist="exponential")
s.mod2<-summary(mod2)

carol_cipro$pred<-predict(mod2)

ggplot(carol_cipro,aes(y=med_cipro,x=ano1,color=sorovar))+
   facet_grid(~material)+
   geom_jitter()+
   geom_line(data = carol_cipro, aes(x = ano1, y = pred))


###############
# Ceftazidime #
###############

carol_ceftaz<-subset(carol,id != 130 & id != 296)


mod3<-survreg(Surv(time=min_ceftaz,time2=max_ceftaz,type="interval2")~ano:sorovar:material+sorovar+material,
              cluster=projeto,data=carol_ceftaz,dist="exponential")
s.mod3<-summary(mod3)

carol_ceftaz$pred<-predict(mod3)

ggplot(carol_ceftaz,aes(y=med_ceftaz,x=ano1,color=sorovar))+
   facet_grid(~material)+
   geom_jitter()+
   geom_line(data = carol_ceftaz, aes(x = ano1, y = pred))

#############
# Ceftaxime #
#############

carol_cefo<-subset(carol,id != 113 & id!=250 & id!=82 & id != 296)

mod4<-survreg(Surv(time=min_cefo,time2=max_cefo,type="interval2")~ano:sorovar:material+sorovar+material,
              cluster=projeto,data=carol_cefo,dist="exponential")
s.mod4<-summary(mod4)


carol_cefo$pred<-predict(mod4)

ggplot(carol_cefo,aes(y=med_cefo,x=ano1,color=sorovar))+
   facet_grid(~material)+
   geom_jitter()+
   geom_line(data = carol_cefo, aes(x = ano1, y = pred))



## All plots together

#Creating a label for ATB and crrecting names of variables

carol_cefo$Atb<-rep("Cefotaxime",dim(carol_cefo)[1])
carol_cefo<-carol_cefo[,c(1:6,25,27,28)]
colnames(carol_cefo)[7]<-"med"

carol_ceftaz$Atb<-rep("Ceftazidime",dim(carol_ceftaz)[1])
carol_ceftaz<-carol_ceftaz[,c(1:6,21,27,28)]
colnames(carol_ceftaz)[7]<-"med"

carol_cipro$Atb<-rep("Ciprofloxacin",dim(carol_cipro)[1])
carol_cipro<-carol_cipro[,c(1:6,17,27,28)]
colnames(carol_cipro)[7]<-"med"


carol_coli$Atb<-rep("Colistin",dim(carol_coli)[1])
carol_coli<-carol_coli[,c(1:6,13,27,28)]
colnames(carol_coli)[7]<-"med"


plot_geral<-rbind.data.frame(carol_cefo,
                             carol_ceftaz,
                             carol_cipro
                             ,carol_coli)

plot_geral$sorovar<-ifelse(plot_geral$sorovar=="Typh","Typhimurium",plot_geral$sorovar)

plot_geral$Atb<- factor(plot_geral$Atb,      # Reordering group factor levels
                         levels = c("Ciprofloxacin", "Colistin", "Ceftazidime", "Cefotaxime"))

mat.labs <- c("Animal", "Food")
names(mat.labs) <- c("Carc", "Food")

ggplot(plot_geral,aes(x=ano,y=med,col=sorovar))+
   theme_minimal()+
   theme(text=element_text(size=20))+
   facet_grid(Atb~material,scales = "free_y",labeller=labeller(material=mat.labs) )+
   geom_jitter()+
   geom_line(data = plot_geral, aes(x = ano, y = pred),size=1)+
   scale_x_continuous(breaks = seq(0,15,2))+
   labs(x ="Year", y = "MIC (mg/L)")+
   scale_color_discrete(name = "Serovar",
                      labels = c("Derby", "Typhimurium"))+
   theme(panel.spacing.x = unit(1, "lines"))+
   theme(panel.spacing.y = unit(1.5, "lines"))
ggsave(here("Figures",'MIC.png'), dpi = 300, height = 10, width = 15, unit = 'in')


#Export tables

final_res<-list(
coli=round(cbind(s.mod1$table[,c(1,5)],confint(mod1)),2),
cipro=round(cbind(s.mod2$table[,c(1,5)],confint(mod2)),2),
ceftaz=round(cbind(s.mod3$table[,c(1,5)],confint(mod3)),2),
cefo=round(cbind(s.mod4$table[,c(1,5)],confint(mod4)),2)
)


sink(here("Outputs","MIC_reg.txt"))
print(final_res)
sink() 
