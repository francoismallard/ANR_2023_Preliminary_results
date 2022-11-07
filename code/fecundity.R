rm(list=ls())
library(emmeans)
all_fec <- list.files("data/fec")

data_fec=NULL
for(i in 1:length(all_fec)){
temp = read.csv(paste0("data/fec/",all_fec[i]),sep=";",h=TRUE,row.names=1)
if(names(temp)[1]=="A" & substring(all_fec[i],1,10)!='2022-10-03'){ # 1rst experiment
  data_fec= rbind(data_fec,data.frame(date=substring(all_fec[i],1,10),
             population=rep(row.names(temp),each=6),
             replicate=rep(colnames(temp[,1:6]),11),
             fec = c(as.matrix(t(temp[,1:6]),1,66))))  
}else if(substring(all_fec[i],1,10)=='2022-10-03'){
  temp1 = temp[,1:6]
  data_fec= rbind(data_fec,data.frame(date=substring(all_fec[i],1,10),
                                      population=rep(row.names(temp1),each=6),
                                      replicate=rep(colnames(temp1[,1:6]),11),
                                      fec = c(as.matrix(t(temp1[,1:6]),1,66))))
  temp1 = temp[,7:14]
  data_fec= rbind(data_fec,data.frame(date=substring(all_fec[i],1,10),
                                      population=rep(row.names(temp1),each=8),
                                      replicate=rep(colnames(temp1[,1:8]),11),
                                      fec = c(as.matrix(t(temp1[,1:8]),1,88))))
  
}else{
  data_fec= rbind(data_fec,data.frame(date=substring(all_fec[i],1,10),
                                      population=rep(row.names(temp),each=8),
                                      replicate=rep(colnames(temp[,1:8]),11),
                                      fec = c(as.matrix(t(temp[,1:8]),1,88))))
  
  
}
}

data_fec$evolution = "ANC"
data_fec$evolution[data_fec$population%in%paste0("AA",1:8)] = "AAG"

survival <- read.table('output_files/txt/survival_all_fixed.txt',h=TRUE,sep='\t')
survival=subset(survival,age>2)
dim(data_fec)
dim(merge(data_fec,survival))

### duration interval
duration <- data.frame(date=unique(data_fec$date),duration=c(1,1,1,2,1,1,1,1,1))

data_fec=merge(data_fec,survival)
data_fec=merge(data_fec,duration)
data_fec$fec2 = data_fec$fec/data_fec$alive_herm/data_fec$duration

data_fec$block="B2"
data_fec$block[data_fec$replicate%in%c("A","B","C","D","E","F")] = "B1"

data_fec$treatment="herm"
data_fec$treatment[data_fec$replicate%in%c("D","E","F","K","L","M","N")] = "mixed"
data_fec$treatment[data_fec$PR=="AA4_K"] = "herm"

boxplot(data_fec$fec2~data_fec$evolution)

temp = glm(fec ~ offset(alive_herm) + offset(duration) +evolution,family='quasipoisson',data=data_fec)
#temp2 = glm(data_fec$fec ~  data_fec$evolution,family='quasipoisson')
em_area <- emmeans(temp,~evolution)
#em_area2 <- emmeans(temp2,~evolution)
pairs(em_area)
#pairs(em_area2)

#### Second Pool
temp = glm(fec ~ offset(alive_herm) + offset(duration) + evolution*date-1, family='quasipoisson',data=subset(data_fec,replicate%in%c("G","H","I","J")))
summary(temp)
em_area <- emmeans(temp,~evolution)
pairs(em_area)

temp_2s <- subset(data_fec,replicate%in%c("G","H","I","J"))
new_pred <- unique(temp_2s[,c("evolution","date","age","duration")])
new_pred$alive_herm=10
new_pred = cbind(new_pred,predict(temp,new_pred,se.fit=TRUE,type="response")$fit,predict(temp,new_pred,se.fit=TRUE,type="response")$se.fit)
names(new_pred)[6:7] = c("fit","se.fit")

vcol <- c(rgb(0,0,0,alpha = .7),rgb(.8,.15,.15, alpha = .7))
plot(new_pred$age,new_pred$fit/10,col=c("black",'firebrick3'),pch=16,las=1,bty="n",ylab="Mean fecundity (/day/herm)",ylim=c(0,80),type="n",xlab="Age (days)")
arrows(new_pred$age,(new_pred$fit+new_pred$se.fit)/10,new_pred$age,(new_pred$fit-new_pred$se.fit)/10,angle=90,code=3,length=.1,col=vcol)
points(new_pred$age,new_pred$fit/10,col=c("black",'firebrick3'),pch=16)


#### First Pool - HERM
temp = glm(fec ~ offset(alive_herm) + offset(duration) +  evolution*date-1, family='quasipoisson',data=subset(data_fec,replicate%in%c("A","B","C") & age<9))
summary(temp)
em_area <- emmeans(temp,~evolution)
pairs(em_area)

temp_2s <- subset(data_fec,replicate%in%c("A","B","C") & age<9)
new_pred <- unique(temp_2s[,c("evolution","date","age","duration")])
new_pred$alive_herm=10
new_pred = cbind(new_pred,predict(temp,new_pred,se.fit=TRUE,type="response")$fit,predict(temp,new_pred,se.fit=TRUE,type="response")$se.fit)
names(new_pred)[6:7] = c("fit","se.fit")

vcol <- c(rgb(0,0,0,alpha = .7),rgb(.8,.15,.15, alpha = .7))
plot(new_pred$age-1,new_pred$fit/10,col=c("black",'firebrick3'),pch=16,las=1,bty="n",ylab="Mean fecundity (/day/herm)",type="n",xlab="Age (days)",ylim=c(0,1000))
arrows(new_pred$age-1,(new_pred$fit+new_pred$se.fit)/10,new_pred$age-1,(new_pred$fit-new_pred$se.fit)/10,angle=90,code=3,length=.1,col=vcol)
points(new_pred$age-1,new_pred$fit/10,col=c("black",'firebrick3'),pch=16)


#### First Pool - MIXED
temp = glm(fec ~ offset(alive_herm) + offset(duration) +  evolution*date, family='quasipoisson',data=subset(data_fec,replicate%in%c("D","E","F") & age<8))
summary(temp)
em_area <- emmeans(temp,~evolution*date)
pairs(em_area)

temp_2s <- subset(data_fec,replicate%in%c("D","E","F") & age<8)
new_pred <- unique(temp_2s[,c("evolution","date","age","duration")])
new_pred$alive_herm=5
new_pred = cbind(new_pred,predict(temp,new_pred,se.fit=TRUE,type="response")$fit,predict(temp,new_pred,se.fit=TRUE,type="response")$se.fit)
names(new_pred)[6:7] = c("fit","se.fit")

vcol <- c(rgb(0,0,0,alpha = .7),rgb(.8,.15,.15, alpha = .7))
#plot(new_pred$age-1,new_pred$fit/5,col=c("black",'firebrick3'),pch=16,las=1,bty="n",ylab="Mean fecundity (/day/herm)",type="n",xlab="Age (days)",ylim=c(0,400))
arrows(new_pred$age-.8,(new_pred$fit+new_pred$se.fit)/5,new_pred$age-.8,(new_pred$fit-new_pred$se.fit)/5,angle=90,code=3,length=.1,col=vcol)
points(new_pred$age-.8,new_pred$fit/5,col=c("black",'firebrick3'),pch=8)

test= subset(data_fec,block=="B1" & age<8)
boxplot(test$fec2~test$treatment+test$date,col=c("black","red"))

pdf("plots/Fecundity_summary.pdf",h=12,w=5)
par(mfrow=c(2,1))
#### First Pool - ALL
temp = glm(fec2 ~  evolution*date*treatment,data=subset(data_fec,block=='B1' & age<8))
summary(temp)
em_area <- emmeans(temp,~evolution*treatment)
pairs(em_area)

temp_2s <- subset(data_fec,block=='B1' & age<8)
new_pred <- unique(temp_2s[,c("evolution","date","age","duration","treatment")])
new_pred = cbind(new_pred,predict(temp,new_pred,se.fit=TRUE,type="response")$fit,predict(temp,new_pred,se.fit=TRUE,type="response")$se.fit)
names(new_pred)[6:7] = c("fit","se.fit")




vcol <- c(rgb(0,0,0,alpha = .7),rgb(.8,.15,.15, alpha = .7))
plot(new_pred$age-1,new_pred$fit,col=c("black",'firebrick3'),pch=16,las=1,bty="n",ylab="Mean fecundity (/day/herm)",type="n",xlab="Age (days)",ylim=c(0,400),xlim=c(2.5,7.5))
arrows(new_pred$age - c(1,.8),(new_pred$fit+new_pred$se.fit),new_pred$age- c(1,.8),(new_pred$fit-new_pred$se.fit),angle=90,code=3,length=.1,col=rep(vcol,each=2))
points(new_pred$age - c(1,.8),new_pred$fit,col=rep(c("black",'firebrick3'),each=2),pch=c(16,8))

#### Second Pool - ALL
temp = glm(fec2 ~  evolution*date*treatment,data=subset(data_fec,block=='B2'))
summary(temp)
em_area <- emmeans(temp,~evolution*treatment)
pairs(em_area)

temp_2s <- subset(data_fec,block=='B2')
new_pred <- unique(temp_2s[,c("evolution","date","age","duration","treatment")])
new_pred = cbind(new_pred,predict(temp,new_pred,se.fit=TRUE,type="response")$fit,predict(temp,new_pred,se.fit=TRUE,type="response")$se.fit)
names(new_pred)[6:7] = c("fit","se.fit")

vcol <- c(rgb(0,0,0,alpha = .7),rgb(.8,.15,.15, alpha = .7))
plot(new_pred$age-1,new_pred$fit/5,col=c("black",'firebrick3'),pch=16,las=1,bty="n",ylab="Mean fecundity (/day/herm)",type="n",xlab="Age (days)",ylim=c(0,180),xlim=c(2.5,7.5))
arrows(new_pred$age + c(0,.2),(new_pred$fit+new_pred$se.fit),new_pred$age + c(0,.2),(new_pred$fit-new_pred$se.fit),angle=90,code=3,length=.1,col=rep(vcol,each=2))
points(new_pred$age + c(0,.2),new_pred$fit,col=rep(c("black",'firebrick3'),each=2),pch=c(16,8))
dev.off()


#data_fec$pop2 = data_fec$population
#data_fec$pop2[data_fec$pop2$evolution=='ANC'] = "A00"

#pdf("plots/Fecundity_summary_per_Replicate.pdf",h=12,w=5)
par(mfrow=c(2,1))
#### First Pool - ALL
temp = glm(fec2 ~  population*date*treatment,data=subset(data_fec,block=='B1' & age<8))
summary(temp)
em_area <- emmeans(temp,~evolution*treatment)
pairs(em_area)

temp_2s <- subset(data_fec,block=='B1' & age<8)
new_pred <- unique(temp_2s[,c("population","date","age","duration","treatment")])
new_pred = cbind(new_pred,predict(temp,new_pred,se.fit=TRUE,type="response")$fit,predict(temp,new_pred,se.fit=TRUE,type="response")$se.fit)
names(new_pred)[6:7] = c("fit","se.fit")

vcol <- c(rgb(0,0,0,alpha = .7),rgb(.8,.15,.15, alpha = .7))
v_xpos = rep(c(-.1,0,.1,rep(.5,8)),each=2) + rep(c(-.1,.1),11)+.6
plot(new_pred$age-1,new_pred$fit,col=c("black",'firebrick3'),pch=16,las=1,bty="n",ylab="Mean fecundity (/day/herm)",type="n",xlab="Age (days)",ylim=c(0,400),xlim=c(2.5,7.5))
arrows(new_pred$age - v_xpos,(new_pred$fit+new_pred$se.fit),new_pred$age- v_xpos,(new_pred$fit-new_pred$se.fit),angle=90,code=3,length=.1,col = rep(vcol,c(6,16)))
points(new_pred$age - v_xpos,new_pred$fit,col=rep(c("black",'firebrick3'),c(6,16)),pch=c(16,8))

#### Second Pool - ALL
temp = glm(fec2 ~  population*date*treatment,data=subset(data_fec,block=='B2'))

temp_2s <- subset(data_fec,block=='B2')
new_pred <- unique(temp_2s[,c("population","date","age","duration","treatment")])
new_pred = cbind(new_pred,predict(temp,new_pred,se.fit=TRUE,type="response")$fit,predict(temp,new_pred,se.fit=TRUE,type="response")$se.fit)
names(new_pred)[6:7] = c("fit","se.fit")


vcol <- c(rgb(0,0,0,alpha = .7),rgb(.8,.15,.15, alpha = .7))
v_xpos = rep(c(-.1,0,.1,rep(.5,8)),each=2) + rep(c(-.1,.1),11)
plot(new_pred$age,new_pred$fit,col=c("black",'firebrick3'),pch=16,las=1,bty="n",ylab="Mean fecundity (/day/herm)",type="n",xlab="Age (days)",ylim=c(0,180),xlim=c(2,8))
arrows(new_pred$age - v_xpos,(new_pred$fit+new_pred$se.fit),new_pred$age - v_xpos,(new_pred$fit-new_pred$se.fit),angle=90,code=3,length=.1,col = rep(vcol,c(6,16)))
new_pred$xpos=v_xpos

for(i in unique(new_pred$population)){
  for(j in unique(new_pred$treatment)){

    tempX=subset(new_pred,population==i & treatment==j)
    if(!unique(tempX$population)%in%c("A00r1","A00r2","A00r3"))  points(tempX$age - tempX$xpos,tempX$fit,type="b",col=ifelse(j=='herm',"blue","black"),pch=16)
  }
}

#dev.off()


#### PLOT FOR ANR


library(lme4)
#### First Pool - ALL
mod_P1 = lmer(fec2 ~  evolution:date:treatment-1+(1|population),data=subset(data_fec,block=='B1' & age<8))
temp_P1 <- subset(data_fec,block=='B1' & age<8)
new_pred <- unique(temp_P1[,c("evolution","date","age","duration","treatment")])
new_pred=new_pred[order(new_pred$treatment,new_pred$date,new_pred$evolution),]
confint_P1 = confint(mod_P1)
pred_values_P1 <- cbind(new_pred,predict(mod_P1,new_pred,re.form=NA,se.fit=TRUE))
pred_values_P1=cbind(pred_values_P1,confint_P1[3:nrow(confint_P1),])
names(pred_values_P1)[6:8] = c("fit","low","high")
pred_values_P1$low[pred_values_P1$low<0]=0

#### Second Pool - ALL
mod_P2 = lmer(fec2 ~  evolution:date:treatment-1+(1|population),data=subset(data_fec,block=='B2'))
temp_P2 <- subset(data_fec,block=='B2')
new_pred <- unique(temp_P2[,c("evolution","date","age","duration","treatment")])
new_pred=new_pred[order(new_pred$treatment,new_pred$date,new_pred$evolution),]
confint_P2 = confint(mod_P2)
pred_values_P2 <- cbind(new_pred,predict(mod_P2,new_pred,re.form=NA,se.fit=TRUE))
pred_values_P2=cbind(pred_values_P2,confint_P2[3:nrow(confint_P2),])
names(pred_values_P2)[6:8] = c("fit","low","high")
pred_values_P2$low[pred_values_P2$low<0]=0

vcol_P1 <- c(rgb(0,0,0,alpha = .7),rgb(.8,.15,.15, alpha = .7))[2:1]
vcol_P2 <- c(rgb(0,0,0,alpha = .7),rgb(.8,.15,.15, alpha = .7))[2:1]

pdf("plots/Fec_ANR.pdf",h=5,w=5)
x_pos_P1 = - rep(c(1.25,.75),each=2) + c(-.1,0.1)[as.numeric(as.factor(pred_values_P1$evolution))]
x_pos_P2 =  rep(c(4.23,3.75),each=2) + c(-.1,0.1)[as.numeric(as.factor(pred_values_P1$evolution))]
plot(pred_values_P1$age-1,pred_values_P1$fit,col=c("black",'firebrick3'),pch=16,las=1,bty="n",ylab="Mean fecundity (/day/ind)",type="n",xlab="Age (days)",ylim=c(0,380),xlim=c(2.5,10),xaxt="n")
axis(side=1,at=c(3,4,5,7,8,9),labels=c(3,4,5,3,4,5))
arrows(pred_values_P1$age + x_pos_P1,(pred_values_P1$low),pred_values_P1$age + x_pos_P1,(pred_values_P1$high),angle=90,code=3,length=.05,col=vcol_P1[as.numeric(as.factor(pred_values_P1$evolution))])
points(pred_values_P1$age + x_pos_P1,pred_values_P1$fit,col=c('firebrick3',"black")[as.numeric(as.factor(pred_values_P1$evolution))],pch=c(16,21)[as.numeric(as.factor(pred_values_P1$treatment))])

pred_values_P2 = subset(pred_values_P2,age<=5)
arrows(pred_values_P2$age + x_pos_P2,(pred_values_P2$low),pred_values_P2$age + x_pos_P2,(pred_values_P2$high),angle=90,code=3,length=.05,col=vcol_P2[as.numeric(as.factor(pred_values_P2$evolution))])
points(pred_values_P2$age + x_pos_P2,pred_values_P2$fit,col=c('firebrick3',"black")[as.numeric(as.factor(pred_values_P2$evolution))],pch=c(16,21)[as.numeric(as.factor(pred_values_P2$treatment))])

legend(6,400,bty="n",c("Herm. alone","Herm. + Males"),pch=c(16,21),ncol=1)
legend(6,300,bty="n",c("Ancestral","Evolved"),lwd=1,col=c("black","red"),ncol=1)

dev.off()



