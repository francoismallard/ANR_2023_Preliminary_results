rm(list=ls())
library(emmeans)
all_surv <- list.files("data/surv/")

data_surv=NULL
for(i in 1:length(all_surv)){
  temp= read.csv(paste0("data/surv/",all_surv[i]),sep=";",h=TRUE)
  temp$date=as.Date(substring(all_surv[i],1,10))
  data_surv=rbind(data_surv,temp)
}


data_surv$evolution = "ANC"
data_surv$evolution[data_surv$population%in%paste0("AA",1:8)] = "AAG"
data_surv$PR=paste(data_surv$population,data_surv$replicate,sep="_")
data_surv=subset(data_surv,PR!="_")
data_surv$PR=as.factor(as.character(data_surv$PR))

temp = tapply(data_surv$date,data_surv$PR,min)
temp=data.frame(PR=names(temp),min_date=as.numeric(temp))
data_surv = merge(data_surv,temp)
rm(temp)
data_surv$age = as.numeric(data_surv$date) - data_surv$min_date + 3
table(data_surv$age)

data_surv$max_herm = 10; data_surv$max_male = 0
data_surv$max_male[data_surv$replicate%in%c("D","E","F","K","L","M","N")] = 5
data_surv$max_herm[data_surv$replicate%in%c("D","E","F","K","L","M","N")] = 5
## One exception
data_surv$max_male[data_surv$PR=='AA4_K'] = 0
data_surv$max_herm[data_surv$PR=='AA4_K'] = 10

### Introduce one extra line at age 0-2 for each PR
for(i in unique(data_surv$PR)){
  temp=subset(data_surv,PR==i & age==3)
  temp$age=2;temp$alive_herm=temp$max_herm;temp$alive_male=temp$max_male
  data_surv=rbind(data_surv,temp)
  #temp$age=0;temp$alive_herm=temp$max_herm;temp$alive_male=temp$max_male
  #data_surv=rbind(data_surv,temp)
  
  }
data_surv=data_surv[order(data_surv$PR,data_surv$age),]

### save the table for the fecundity analysis
#write.table(data_surv,file='output_files/txt/survival_all.txt',quote=FALSE,row.names=FALSE,sep="\t")

data_surv$block = "B1" 
data_surv$block[data_surv$replicate%in%c("G","H","I","J","K","L","M","N")] = "B2" 
data_surv$treatment = "herm"
data_surv$treatment[data_surv$replicate%in%c("D","E","F","K","L","M","N")] = "mixed" 
data_surv$treatment[data_surv$PR=="AA4_K"] = "herm"


temp1=subset(data_surv,PR=="AA1_A")
plot(temp1$age,temp1$alive_herm,type="n")
for(i in unique(data_surv$PR)){
  temp1=subset(data_surv,PR==i)
  lines(temp1$age,temp1$alive_herm,type="l")
}

### First batch - pool
d1=subset(data_surv,replicate%in%c("A","B","C"))
temp1=subset(d1,PR=="AA1_A")
plot(temp1$age,temp1$alive_herm,type="n",ylim=c(0,1),las=1,bty="n")
for(i in unique(d1$population)){
  d2=subset(d1,population==i)
  temp= tapply(d2$alive_herm/d2$max_herm,d2$age,mean)
  vcol = ifelse(substring(i,1,2)=="AA","red","black")
  lines(as.numeric(names(temp)),as.numeric(temp),type="l",col=vcol)
}

#d1=subset(data_surv,replicate%in%c("D","E","F"))
#for(i in unique(d1$population)){
#  d3=subset(d1,population==i)
#  d3$alive_herm[is.na(d3$alive_herm)] = 0
#  temp= tapply(d3$alive_herm/d3$max_herm,d3$age,mean)
#  vcol = ifelse(substring(i,1,2)=="AA","yellow","green")
#  lines(as.numeric(names(temp)),as.numeric(temp),type="l",col=vcol)
#}

### Second batch - pool
d1=subset(data_surv,replicate%in%c("G","H","I","J") | PR=='AA4_K')
d2=subset(data_surv,!(replicate%in%c("G","H","I","J") | PR=='AA4_K') & block=='B2')

temp1=subset(d1,PR=="AA1_G")
plot(temp1$age,temp1$alive_herm,type="n",ylim=c(0,1),las=1,bty="n",xlab="Age (days)",ylab="Proportion alive (n=80)")
for(i in unique(d1$population)){
  d3=subset(d1,population==i)
  d3$alive_herm[is.na(d3$alive_herm)] = 0
  temp= tapply(d3$alive_herm/d3$max_herm,d3$age,mean)
  #temp= tapply(d3$alive_herm,d3$age,sum)
  vcol = ifelse(substring(i,1,2)=="AA","red","black")
  lines(as.numeric(names(temp)),as.numeric(temp),type="l",col=vcol)
}

for(i in unique(d2$population)){
  d3=subset(d2,population==i)
  d3$alive_herm[is.na(d3$alive_herm)] = 0
  temp= tapply(d3$alive_herm/d3$max_herm,d3$age,mean)
  vcol = ifelse(substring(i,1,2)=="AA","yellow","green")
  lines(as.numeric(names(temp)),as.numeric(temp),type="l",col=vcol)
}


#######


data_surv$nb_morts_herm = 0

data_surv$nb_morts_male = 0
data_surv2=NULL
for(i in unique(data_surv$PR)){
  temp=subset(data_surv, PR==i)
  temp$nb_morts_herm[2:nrow(temp)] = temp$alive_herm[1:I(nrow(temp)-1)] - temp$alive_herm[2:nrow(temp)]
  temp$nb_morts_male[2:nrow(temp)] = temp$alive_male[1:I(nrow(temp)-1)] - temp$alive_male[2:nrow(temp)]
  data_surv2=rbind(data_surv2,temp)
}
data_surv=data_surv2
rm(data_surv2);gc()

#### !!!!!
subset(data_surv,nb_morts_herm<0)
data_surv$nb_morts_herm[data_surv$nb_morts_herm<0] = 0

library(survival)

data_for_cox <- NULL
for(i in unique(data_surv$PR)){
  
  temp=subset(data_surv, PR==i)
  temp$nb_morts_herm[is.na(temp$nb_morts_herm)]=0
  if(sum(temp$nb_morts_herm)>0){
  temp_df=data.frame(age=rep(temp$age,temp$nb_morts_herm),status=2)
  }
  #if(sum(temp$nb_disparus>0)){
  #  temp_df=rbind(temp_df,data.frame(age=rep(temp$age,temp$nb_disparus),status=1))
  #}
  temp_df$PR=i	
  data_for_cox=rbind(data_for_cox, temp_df)
}

head(data_surv)
data_for_cox = merge(data_for_cox,unique(data_surv[,c("PR","population","block","replicate","evolution","treatment")]))
head(data_for_cox)

str(data_for_cox)
res.cox <- coxph(Surv(age, status) ~ evolution, data = subset(data_for_cox,treatment=="herm" & block=="B1"))
summary(res.cox)

res.cox <- coxph(Surv(age, status) ~ evolution, data = subset(data_for_cox,treatment=="herm" & block=="B2"))
summary(res.cox)
plot(survfit(Surv(age, status) ~ evolution, data=subset(data_for_cox,treatment=="herm" & block=="B1")))

res.cox <- coxph(Surv(age, status) ~ evolution*treatment, data = subset(data_for_cox, block=="B2"))
summary(res.cox)


plot(survfit(Surv(age, status) ~ evolution, data=subset(data_for_cox,treatment=="herm" & block=="B1")),xlim=c(0,20))
par(new=FALSE)
lines(survfit(Surv(age, status) ~ evolution, data=subset(data_for_cox,treatment=="herm" & block=="B2")),xlim=c(0,20),col='red')

library(ggfortify)
autoplot(survfit(Surv(age, status) ~ evolution+block, data=subset(data_for_cox,treatment=="herm")))

confint(survfit(Surv(age, status) ~ evolution, data=subset(data_for_cox,treatment=="herm" & block=="B2")))
plot(res.cox)

plot(survfit(Surv(age, status)~ population, subset(data_for_cox,treatment=="herm" & block=="B1")),col=rep(c("black","red"),c(3,8)))

data_P1 <- subset(data_for_cox,treatment=="herm" & block=="B1")
t1=tapply(data_P1$age,data_P1$evolution,mean)
t1_sd=tapply(data_P1$age,data_P1$evolution,sd)

data_P2 <- subset(data_for_cox,treatment=="herm" & block=="B2")
t2=tapply(data_P2$age,data_P2$evolution,mean)
t2_sd=tapply(data_P2$age,data_P2$evolution,sd)

barplot(c(t1,t2),col=c("black","firebrick3","black","firebrick3"),space=c(0,0,1,0),ylab="Mean age at death",ylim=c(0,20))
arrows(c(.5,1.5,3.5,4.5),c(t1-t1_sd,t2-t2_sd),c(.5,1.5,3.5,4.5),c(t1+t1_sd,t2+t2_sd),
       code=3,length=.05,angle=90)



res.cox <- coxph(Surv(age, status) ~ evolution , data = subset(data_for_cox,treatment=="herm" & block=="B1"))
autoplot(survfit(Surv(age, status) ~ evolution, data=subset(data_for_cox,treatment=="herm" & block=="B2")))
summary(res.cox)
res.cox <- coxph(Surv(age, status) ~ evolution, data = subset(data_for_cox,treatment=="herm" & block=="B2"))
summary(res.cox)



mod_surv <- lmer(age~evolution+(1|population),data=subset(data_for_cox,treatment=="herm" & block=="B1"))
mod_surv2 <- lmer(age~1+(1|population),data=subset(data_for_cox,treatment=="herm" & block=="B1"))
anova(mod_surv,mod_surv2)


mod_surv <- glmer(age~ evolution+(1|population),data=subset(data_for_cox,treatment=="herm" & block=="B2"),family="poisson")
mod_surv2 <- glmer(age~ 1+(1|population),data=subset(data_for_cox,treatment=="herm" & block=="B2"),family="poisson")
anova(mod_surv,mod_surv2)


summary(mod_surv)


### First batch - pool
pdf("plots/Surv_ANR.pdf",h=4,w=6)
par(mfrow=c(1,2))
d1=subset(data_surv,replicate%in%c("A","B","C"))
temp1=subset(d1,PR=="AA1_A")
plot(temp1$age,temp1$alive_herm,type="n",ylim=c(0,1),las=1,bty="n",xlab="Age (days)",ylab="Proportion alive (n=80)")
for(i in unique(d1$population)){
  d2=subset(d1,population==i)
  temp= tapply(d2$alive_herm/d2$max_herm,d2$age,mean)
  vcol = ifelse(substring(i,1,2)=="AA","red","black")
  lines(as.numeric(names(temp)),as.numeric(temp),type="l",col=vcol)
}
### Second batch - pool
d1=subset(data_surv,replicate%in%c("G","H","I","J") | PR=='AA4_K')
d2=subset(data_surv,!(replicate%in%c("G","H","I","J") | PR=='AA4_K') & block=='B2')

temp1=subset(d1,PR=="AA1_G")
plot(temp1$age,temp1$alive_herm,type="n",ylim=c(0,1),las=1,bty="n",xlab="Age (days)",ylab="Proportion alive (n=80)")
for(i in unique(d1$population)){
  d3=subset(d1,population==i)
  d3$alive_herm[is.na(d3$alive_herm)] = 0
  temp= tapply(d3$alive_herm/d3$max_herm,d3$age,mean)
  #temp= tapply(d3$alive_herm,d3$age,sum)
  vcol = ifelse(substring(i,1,2)=="AA","red","black")
  lines(as.numeric(names(temp)),as.numeric(temp),type="l",col=vcol)
}
legend(9,1,c("Ancestral","Evolved"),lwd=c(1,1),col=c("black","red"),bty="n",cex=.8)
dev.off()





