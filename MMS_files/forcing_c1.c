  /*
  Version: 1.1
  */
  t2 = 3.141592653589793*3.141592653589793;
  t3 = epsilon*epsilon;
  t4 = 1.0/Lx;
  t6 = 1.0/Ly;
  t8 = 1.0/Lz;
  t10 = 1.0/T;
  t11 = 1.0/kb;
  t12 = 1.0/na;
  t13 = 1.0/tort;
  t5 = t4*t4;
  t7 = t6*t6;
  t9 = t8*t8;
  t14 = t6*y*3.141592653589793;
  t16 = t4*x*3.141592653589793*2.0;
  t17 = t8*z*3.141592653589793*2.0;
  t18 = t*t10*3.141592653589793*2.0;
  t15 = cos(t14);
  t19 = cos(t16);
  t20 = cos(t17);
  t21 = sin(t16);
  t22 = sin(t18);
  t23 = t15+1.0;
  t24 = t20*t20;
  t25 = t22*t22;
  t26 = t23*t23;
  t27 = (c2Eq*epsilon*t20*t21*t22*t23)/2.0;
  t28 = c1Eq+t27;
  t0 = coeff_hrxn*poro*(hRxnRate7*t28-hRxnRate6*(c5Eq+(c5Eq*epsilon*t19*t20*t22*t23)/2.0)-hRxnRate8*(c4Eq+(c4Eq*epsilon*t19*t20*t22*t23)/2.0)+hRxnRate5*t28*(c7Eq+(c7Eq*epsilon*t19*t20*t22*t23)/2.0))+D1*coeff_diff*poro*t13*((c2Eq*epsilon*t2*t7*t15*t20*t21*t22)/2.0+c2Eq*epsilon*t2*t5*t20*t21*t22*t23*2.0+c2Eq*epsilon*t2*t9*t20*t21*t22*t23*2.0)+D1*coeff_migr*e*poro*t10*t11*t13*z1*((epsilon*phiEq*t2*t7*t15*t19*t20*t22*t28)/2.0+epsilon*phiEq*t2*t5*t19*t20*t22*t23*t28*2.0+epsilon*phiEq*t2*t9*t19*t20*t22*t23*t28*2.0-(c2Eq*phiEq*t2*t3*t7*t19*t21*t24*t25*pow(sin(t14),2.0))/4.0-c2Eq*phiEq*t2*t3*t9*t19*t21*t25*t26*pow(sin(t17),2.0)+c2Eq*phiEq*t2*t3*t5*t19*t21*t24*t25*t26)-(coeff_frxn*fRxn1MoleReactant*fRxn1ExchangeCurrent*fRxn1SurfaceRoughness*t12*t28*exp(-e*fRxn1TransferCoefficient*t10*t11*(fRxn1StandardElectrodePotential+phiEq-phiElectrode+(epsilon*phiEq*t19*t20*t22*t23)/2.0))*exp(-fRxn1ActivationEnergy*t10*t11*t12))/(e*fRxn1CRef)+c2Eq*coeff_time*epsilon*poro*t10*t20*t21*t23*3.141592653589793*cos(t18);
