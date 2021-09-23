FFNGM<-function(x,t, model=c("Monomolecular", "Logistic", "Gompertz"),k,y, r, h){
  nlc <- nls.control(maxiter = 1000)
  grz1=nls(x~k-(k-y)*exp(-r*t), control=nlc, start=list(k=k ,y=y,r=r))# Monomolecular growth model
  ju1=summary(grz1)
  ek1=coef(grz1)[1]
  ey1=coef(grz1)[2]
  er1=coef(grz1)[3]
  grz2=nls(x~k/(1+(k/y-1)* exp(-r*t)),  control=nlc, start=list(k=k,y=y,r=r))#logistic model
  ju2=summary(grz2)
  ek2=coef(grz2)[1]
  ey2=coef(grz2)[2]
  er2=coef(grz2)[3]
  grz3=nls(x~k*exp(log(y/k)* exp(-r*t)), control = nlc, start=list(k=k,y=y,r=r))#Gompertz model
  ek3=coef(grz3)[1]
  ey3=coef(grz3)[2]
  er3=coef(grz3)[3]
  ju3=summary(grz3)
  t1<-length(x)
  hh1<-vector("numeric")
  for (i in 1: t1) {
    hh1[[i]]=ek1-(ek1-ey1)*exp(-er1*t[i])}
  hh2<-vector("numeric")
  for (i in 1: t1) {
    hh2[[i]]=ek2/(1+(ek2/ey2-1)* exp(-er2*t[i]))}
  hh3<-vector("numeric")
  for (i in 1: t1) {
    hh3[[i]]=ek3*exp(log(ey3/ek3)* exp(-er3*t[i]))}
  MAE_Mono=mean(abs(x - hh1))
  MAPE_Mono=(mean(abs(x - hh1)/x))*100
  MSE_Mono=(mean((x - hh1)^2))
  RMSE_Mono=sqrt(mean((x - hh1)^2))
  MAE_Log=mean(abs(x - hh2))
  MAPE_Log=(mean(abs(x - hh2)/x))*100
  MSE_Log=(mean((x - hh2)^2))
  RMSE_Log=sqrt(mean((x - hh2)^2))
  MAE_Gom=mean(abs(x - hh3))
  MAPE_Gom=(mean(abs(x - hh3)/x))*100
  MSE_Gom=(mean((x - hh3)^2))
  RMSE_Gom=sqrt(mean((x - hh3)^2))
  hh4<-vector("numeric")
  for (i in (t1+1): (t1+h)) {
    hh4[[i]]=ek1-(ek1-ey1)*exp(-er1*i)}
  hh44=hh4[(t1+1):(t1+h)]
  hh5<-vector("numeric")
  for (i in (t1+1): (t1+h)) {
    hh5[[i]]=ek2/(1+(ek2/ey2-1)* exp(-er2*i))}
  hh55=hh5[(t1+1):(t1+h)]
  hh6<-vector("numeric")
  for (i in (t1+1): (t1+h)) {
    hh6[[i]]=ek3*exp(log(ey3/ek3)* exp(-er3*i))}
  hh66=hh6[(t1+1):(t1+h)]
  if (model=="Monomolecular") {return(list(moodelsummary=ju1,fitted.values=hh1, MAE=MAE_Mono,MAPE=MAPE_Mono, MSE=MSE_Mono, RMSE=RMSE_Mono,forecated.values=hh44))}
  if (model=="Logistic") {return(list(modelsummary=ju2, fitted.values=hh2,MAE=MAE_Log,MAPE=MAPE_Log, MSE=MSE_Log, RMSE=RMSE_Log, forecasted.values=hh55))}
  if (model=="Gompertz") {return(list(modelsummary=ju3, fitted.values=hh3, MAE=MAE_Gom, MAPE=MAPE_Gom, MSE=MSE_Gom, RMSE=RMSE_Gom, forecasted.values=hh66))}
  else {return("not a valid model")}
}
