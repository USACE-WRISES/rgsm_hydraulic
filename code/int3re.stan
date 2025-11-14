  
data{
 int Nexperts;
 real lNz [4];
 real lCz [4];
 int Nyears;
 int Nstrata;
 int NR0;
 int NR1;
 int NAVsamps;
 int Nobs_mon1;
 int Nobs_mon0;
 int Nobs_mon01;
 int Nobs_monV;
 int Ntotjul;
 real mon1_effQ [Nobs_mon1,2];//effort,Q,
 real mon0_effQ [Nobs_mon0,2];//effort,Q,
 real mon01_effQ [Nobs_mon01,2];//effort,Q,
 real monV_effQ [Nobs_monV,3];//effort,Q,reweighting
 matrix [NAVsamps,4] mAV;//Cstrata,Q,pool,total  //don't actually need first column
 int mon1[Nobs_mon1,5]; //monitoring data for age1+ during summer - year, julian,Cstrata,habitat,catch
 int mon0[Nobs_mon0,5]; //monitoring data for age0..only June through Sept - year, julian,Cstrata,habitat,catch
 int mon01[Nobs_mon01,5]; //monitoring data for age0 +age1..only June through October - year, julian,Cstrata,habitat,catch
 int monV[Nobs_monV,5]; //monitoring data for vie - year, julian,Cstrata,habitat,catch
 real StrataLen [Nyears,Ntotjul,Nstrata];
 real cum_nd [Nyears,Ntotjul,Nstrata];
 real prop_nd [Nyears,Ntotjul,Nstrata];
 real cum_phiR [Nyears,Ntotjul,Nstrata];
 int R0 [NR0,4];
 int R1 [NR1,4];
 matrix [6,Nstrata] w_um;
 matrix [6,Nstrata] nmons; //average number of winter months for each year/strata 
 real w_m [(Nyears-5),Nstrata,2]; //last dimension is whether released ~2 or ~5 months; 
 matrix [Nyears,Nstrata] X;
 real emove [Nexperts,2];
 vector [8] ewidths [Nexperts,3,2];
 real refQ [8];
 }


parameters{
  real <lower=7,upper=13> lN0 [Nstrata,2]; //initial abundance
  //reproduction parameters
  real <lower=300,upper=1500> a; //per spawner production of age 1
  real <lower=0.5,upper=3> beta_stk;
  real <lower=0.5,upper=3> beta_2;
  vector [Nyears] lbeta_eps [Nstrata];
  real mu_lbeta [Nstrata]; //average log carrying capacity
  real <lower=0> sd_lbeta;
  real <lower=-20,upper=20> B_lbeta;
  // survival parameters
  real <lower=-8, upper=-3> mu_lM0[Nstrata];
  real <lower=-8, upper=-3> mu_lM1[Nstrata];
  real <lower=-8, upper=-3> mu_lMw[Nstrata];
  real <lower=-8, upper=-3> lM0 [Nstrata,Nyears]; //natural mortality rate
  real <lower=-8, upper=-3> lM1 [Nstrata,Nyears];
  real <lower=-8, upper=-3> lMw [Nstrata,(Nyears+1)];
  real <lower=0.001,upper=1> sd_lM;
  real <lower=0,upper=1> irphi;
  //move
  real <lower=0,upper=1> move; //movement out of drying reach
  //sampling stuff
  real <lower=0,upper=1> lA0_perpool;
  real <lower=-10,upper=0> AtQ_perpool;
  real <lower=0.01> sigma_obspool;
  real <lower=6,upper=10> logsl_width [Nstrata];
  real <lower=50,upper=200> bankfull [Nstrata];
  real <lower=1,upper=100> alpha0_int;
  real <lower=1,upper=100> alpha1_int;
  real <lower=1,upper=40> alpha0_max;
  real <lower=1,upper=40> alpha1_max;
  real <lower=0.01> sz;
  real <lower=0.01> rsz;
  real <lower=0,upper=1> p1;
  real <lower=0,upper=1> p0;
  real <lower=0,upper=1> rp1; //effective capture probability during rescue - age 1+
  real <lower=0,upper=1> rp0;
}

transformed parameters{
  real <lower=0> talpha0 [2];
  real <lower=0> talpha1 [2];
  //transformed parameters
  matrix [Nyears,Nstrata] Rf;
  real M0 [Nstrata,Nyears];
  real M1 [Nstrata,Nyears];
  real Mw [Nstrata,(Nyears+1)];
  //sampling
  real sl_width [Nstrata];
  real A0_perpool;
  vector [NAVsamps] pool_area;
  vector [8] pred_width [Nstrata];
  vector [Nobs_mon1] lpmon1;
  vector [Nobs_mon0] lpmon0;  
  vector [Nobs_mon01] lpmon01;
  vector [Nobs_monV] lpmonV;
  vector [NR0] lpR0;
  vector [NR1] lpR1;
  vector [4] tN;
  vector [4] lN;
  real totpool;
  real totrun;
  //abundance related
  real Sp_N [(Nyears+1),Nstrata,3];
  real N0 [Nyears,Nstrata];
  real N1 [(Nyears+1),Nstrata];
  real F_N [Nyears,Nstrata,2];
  real effS [Nyears,Nstrata];
  real Nvie [(Nyears-5),Nstrata];
  //
    A0_perpool=logit(lA0_perpool);
	for (k in 1:Nstrata){
		sl_width[k]=exp(logsl_width[k]);
		}
	for (k in 1:Nstrata){
		Rf[,k]=exp(lbeta_eps[k,]+B_lbeta*X[,k]);
		for (t in 1:Nyears){
			M0[k,t]=exp(lM0[k,t]);
			M1[k,t]=exp(lM1[k,t]);
			}
		for (t in 1:(Nyears+1)){
			Mw[k,t]=exp(lMw[k,t]);
			}}
	for(t in 1:(Nyears-5)){
		for (k in 1:Nstrata){
			Nvie[t,k]=(w_m[t,k,1]*exp(-60*Mw[k,(t+6)])+w_m[t,k,2]*exp(-151*Mw[k,(t+6)]))*irphi;
			}}
	//population dynamics
	for (k in 1:Nstrata){
		Sp_N[1,k,1]=exp(lN0[k,1]); // age 1 abundance mar 31
		Sp_N[1,k,2]=exp(lN0[k,2]); // age 2+ abundance mar 31
		N1[1,k]=sum(Sp_N[1,k,1:2]);
		for (t in 1:6){
			Sp_N[t,k,3]=w_um[t,k]*irphi*exp(-Mw[k,t]*nmons[t,k]*30);
			effS[t,k]=Sp_N[t,k,1]+Sp_N[t,k,2]*beta_2+Sp_N[t,k,3]*beta_stk;
			N0[t,k]=a*(effS[t,k])/(1+a*(effS[t,k])/Rf[t,k]); // age 0 abundance on july 1 - day 92 
			F_N[t,k,1]= N0[t,k]*exp(-M0[k,t]*124)*((1-cum_nd[t,91,k])+(move-1)*(cum_nd[t,215,k]-cum_nd[t,91,k])+(1-move)*rp0*(cum_phiR[t,215,k]-cum_phiR[t,91,k]));
			F_N[t,k,2]=N1[t,k]*exp(-M1[k,t]*215)*(1+(move-1)*cum_nd[t,215,k]+(1-move)*rp1*cum_phiR[t,215,k]);
			Sp_N[(t+1),k,1]=F_N[t,k,1]*exp(-Mw[k,(t+1)]*150);
			Sp_N[(t+1),k,2]=F_N[t,k,2]*exp(-Mw[k,(t+1)]*150);
			N1[(t+1),k]=sum(Sp_N[(t+1),k,1:2]);
			}		
		for(t in 7:Nyears){
			Sp_N[t,k,3]=Nvie[(t-6),k];
			effS[t,k]=Sp_N[t,k,1]+Sp_N[t,k,2]*beta_2+Sp_N[t,k,3]*beta_stk;
			N0[t,k]=a*(effS[t,k])/(1+a*(effS[t,k])/Rf[t,k]); // age 0 abundance on july 1 - day 92 
			F_N[t,k,1]= N0[t,k]*exp(-M0[k,t]*124)*((1-cum_nd[t,91,k])+(move-1)*(cum_nd[t,215,k]-cum_nd[t,91,k])+(1-move)*rp0*(cum_phiR[t,215,k]-cum_phiR[t,91,k]));
			F_N[t,k,2]=N1[t,k]*exp(-M1[k,t]*215)*(1+(move-1)*cum_nd[t,215,k]+(1-move)*rp1*cum_phiR[t,215,k]);
			Sp_N[(t+1),k,1]=F_N[t,k,1]*exp(-Mw[k,(t+1)]*150);
			Sp_N[(t+1),k,2]=F_N[t,k,2]*exp(-Mw[k,(t+1)]*150);
			N1[(t+1),k]=sum(Sp_N[(t+1),k,1:2]);
			}}			
	for (i in 1:NAVsamps){
		pool_area[i]=exp(A0_perpool+AtQ_perpool*mAV[i,2])*mAV[i,4];
		}
	for (j in 1:Nstrata){
		for (k in 1:8){
			pred_width[j,k]=sl_width[j]*refQ[k]/(1+sl_width[j]*refQ[k]/bankfull[j]);
			}}	
	tN[1]=sum(F_N[7,,1])+sum(F_N[7,,2]);
	tN[2]=sum(F_N[8,,1])+sum(F_N[8,,2]);
	tN[3]=sum(F_N[9,,1])+sum(F_N[9,,2]);
	tN[4]=sum(F_N[10,,1])+sum(F_N[10,,2]);
	lN[1]=log(tN[1]);
	lN[2]=log(tN[2]);
	lN[3]=log(tN[3]);
	lN[4]=log(tN[4]);
	talpha1[2]=1;
	for (i in 1:Nobs_mon1){
		totpool=exp(A0_perpool+AtQ_perpool*mon1_effQ[i,2])*200*sl_width[mon1[i,3]]*mon1_effQ[i,2]/(1+sl_width[mon1[i,3]]*mon1_effQ[i,2]/bankfull[mon1[i,3]]);
		totrun=(1-exp(A0_perpool+AtQ_perpool*mon1_effQ[i,2]))*200*sl_width[mon1[i,3]]*mon1_effQ[i,2]/(1+sl_width[mon1[i,3]]*mon1_effQ[i,2]/bankfull[mon1[i,3]]);
		talpha1[1]=alpha1_int*mon1_effQ[i,2]/(1+alpha1_int*mon1_effQ[i,2]/alpha1_max);
		lpmon1[i]=log((p1*mon1_effQ[i,1]*talpha1[mon1[i,4]]*.2*N1[mon1[i,1],mon1[i,3]]*
			exp(-M1[mon1[i,3],mon1[i,1]]*mon1[i,2])*(1+(move-1)*cum_nd[mon1[i,1],mon1[i,2],mon1[i,3]]+(1-move)*rp1*cum_phiR[mon1[i,1],mon1[i,2],mon1[i,3]]))/
			(StrataLen[mon1[i,1],mon1[i,2],mon1[i,3]]*(talpha1[1]*totpool+totrun)));
		}
	talpha0[2]=1;
	for (i in 1:Nobs_mon0){
		totpool=exp(A0_perpool+AtQ_perpool*mon0_effQ[i,2])*200*sl_width[mon0[i,3]]*mon0_effQ[i,2]/(1+sl_width[mon0[i,3]]*mon0_effQ[i,2]/bankfull[mon0[i,3]]);
		totrun=(1-exp(A0_perpool+AtQ_perpool*mon0_effQ[i,2]))*200*sl_width[mon0[i,3]]*mon0_effQ[i,2]/(1+sl_width[mon0[i,3]]*mon0_effQ[i,2]/bankfull[mon0[i,3]]);
		talpha0[1]=alpha0_int*mon0_effQ[i,2]/(1+alpha0_int*mon0_effQ[i,2]/alpha0_max);
		lpmon0[i]=log((p0*mon0_effQ[i,1]*talpha0[mon0[i,4]]*.2*N0[mon0[i,1],mon0[i,3]]*exp(-M0[mon0[i,3],mon0[i,1]]*(mon0[i,2]-91))*
			((1-cum_nd[mon0[i,1],91,mon0[i,3]])+(move-1)*(cum_nd[mon0[i,1],mon0[i,2],mon0[i,3]]-cum_nd[mon0[i,1],91,mon0[i,3]])+
			(1-move)*rp0*(cum_phiR[mon0[i,1],mon0[i,2],mon0[i,3]]-cum_phiR[mon0[i,1],91,mon0[i,3]])))/
			(StrataLen[mon0[i,1],mon0[i,2],mon0[i,3]]*(talpha0[1]*totpool+totrun)));
		}
	for (i in 1:Nobs_mon01){
		totpool=exp(A0_perpool+AtQ_perpool*mon01_effQ[i,2])*200*sl_width[mon01[i,3]]*mon01_effQ[i,2]/(1+sl_width[mon01[i,3]]*mon01_effQ[i,2]/bankfull[mon01[i,3]]);
		totrun=(1-exp(A0_perpool+AtQ_perpool*mon01_effQ[i,2]))*200*sl_width[mon01[i,3]]*mon01_effQ[i,2]/(1+sl_width[mon01[i,3]]*mon01_effQ[i,2]/bankfull[mon01[i,3]]);
		talpha0[1]=alpha0_int*mon01_effQ[i,2]/(1+alpha0_int*mon01_effQ[i,2]/alpha0_max);
		talpha1[1]=alpha1_int*mon01_effQ[i,2]/(1+alpha1_int*mon01_effQ[i,2]/alpha1_max);
		lpmon01[i]=log((p0*mon01_effQ[i,1]*talpha0[mon01[i,4]]*.2*N0[mon01[i,1],mon01[i,3]]*exp(-M0[mon01[i,3],mon01[i,1]]*(mon01[i,2]-91))*
			((1-cum_nd[mon01[i,1],91,mon01[i,3]])+(move-1)*(cum_nd[mon01[i,1],mon01[i,2],mon01[i,3]]-cum_nd[mon01[i,1],91,mon01[i,3]])+
			(1-move)*rp0*(cum_phiR[mon01[i,1],mon01[i,2],mon01[i,3]]-cum_phiR[mon01[i,1],91,mon01[i,3]])))/
			(StrataLen[mon01[i,1],mon01[i,2],mon01[i,3]]*(talpha0[1]*totpool+totrun))+
			(p1*mon01_effQ[i,1]*talpha1[mon01[i,4]]*.2*N1[mon01[i,1],mon01[i,3]]*
			exp(-M1[mon01[i,3],mon01[i,1]]*mon01[i,2])*(1+(move-1)*cum_nd[mon01[i,1],mon01[i,2],mon01[i,3]]+(1-move)*rp1*cum_phiR[mon01[i,1],mon01[i,2],mon01[i,3]]))/
			(StrataLen[mon01[i,1],mon01[i,2],mon01[i,3]]*(talpha1[1]*totpool+totrun)));
		}
	for (i in 1:Nobs_monV){
		totpool=exp(A0_perpool+AtQ_perpool*monV_effQ[i,2])*200*sl_width[monV[i,3]]*monV_effQ[i,2]/(1+sl_width[monV[i,3]]*monV_effQ[i,2]/bankfull[monV[i,3]]);
		totrun=(1-exp(A0_perpool+AtQ_perpool*monV_effQ[i,2]))*200*sl_width[monV[i,3]]*monV_effQ[i,2]/(1+sl_width[monV[i,3]]*monV_effQ[i,2]/bankfull[monV[i,3]]);
		talpha1[1]=alpha1_int*monV_effQ[i,2]/(1+alpha1_int*monV_effQ[i,2]/alpha1_max);
		lpmonV[i]=log(monV_effQ[i,3]*(p1*monV_effQ[i,1]*talpha1[monV[i,4]]*Nvie[(monV[i,1]-6),monV[i,3]]*exp(-M1[monV[i,3],monV[i,1]]*monV[i,2]))/
			(talpha1[1]*totpool+totrun));
		}
//add R0,R1
	for (i in 1:NR0){
		lpR0[i]=log(rp0*prop_nd[R0[i,1],R0[i,2],R0[i,3]]*N0[R0[i,1],R0[i,3]]*exp(-M0[R0[i,3],R0[i,1]]*(R0[i,2]-91))*(1-move));
		}
	for (i in 1:NR1){
		lpR1[i]=log(rp1*prop_nd[R1[i,1],R1[i,2],R1[i,3]]*N1[R1[i,1],R1[i,3]]*exp(-M1[R1[i,3],R1[i,1]]*R1[i,2])*(1-move));
		}
	}
  
model{
  //priors
  a~normal(675,135); //derived by multiplying 3000 from caldwell times .6 fertilization rate and assumed survival of 0.75 and sex ratio of 0.5 and assuming overal cv of 0.2
  beta_stk~normal(1,0.1); //assume similar fecundity for stocked and non-stocked
  beta_2~normal(2,0.1); // roughly based on Caldwell
  for (i in 1:Nexperts){
	emove[i,1]~normal(move,emove[i,2]); //derived from question 5a
	ewidths[i,1,1,]~normal(pred_width[1,],ewidths[i,1,2,]); //derived from question 6a
	ewidths[i,2,1,]~normal(pred_width[2,],ewidths[i,2,2,]); //derived from question 6b
	ewidths[i,3,1,]~normal(pred_width[3,],ewidths[i,3,2,]); //derived from question 6c
	}	
  for (k in 1:Nstrata){
    lbeta_eps[k,]~normal(mu_lbeta[k],sd_lbeta);
	for (t in 1:Nyears){
		lM0[k,t]~normal(mu_lM0[k],sd_lM);
		lM1[k,t]~normal(mu_lM1[k],sd_lM);
		}
	for (t in 1:(Nyears+1)){
		lMw[k,t]~normal(mu_lMw[k],sd_lM);
		}	
	}
  //data
  mAV[,3] ~ normal(pool_area,sigma_obspool);
  mon0[,5] ~ neg_binomial_2_log(lpmon0, sz);
  mon1[,5] ~ neg_binomial_2_log(lpmon1, sz);
  mon01[,5] ~ neg_binomial_2_log(lpmon01, sz);
  monV[,5] ~ neg_binomial_2_log(lpmonV, sz);
  R0[,4] ~ neg_binomial_2_log(lpR0, rsz);
  R1[,4] ~ neg_binomial_2_log(lpR1, rsz);
  lNz ~ normal(lN,lCz);
  }

 
