  /*
  Version: 1.1
  */
  t2 = 3.141592653589793*3.141592653589793;
  t3 = epsilon*epsilon;
  t4 = 1.0/Lx;
  t6 = 1.0/Ly;
  t8 = 1.0/Lz;
  t10 = 1.0/T;
  t11 = 1.0/tort;
  t5 = t4*t4;
  t7 = t6*t6;
  t9 = t8*t8;
  t12 = t6*y*3.141592653589793;
  t14 = t4*x*3.141592653589793*2.0;
  t15 = t8*z*3.141592653589793*2.0;
  t16 = t*t10*3.141592653589793*2.0;
  t13 = cos(t12);
  t17 = cos(t14);
  t18 = cos(t15);
  t19 = sin(t14);
  t20 = sin(t16);
  t21 = t13+1.0;
  t22 = t17*t17;
  t23 = t18*t18;
  t24 = t20*t20;
  t25 = t21*t21;
  t26 = (c3Eq*epsilon*t17*t18*t20*t21)/2.0;
  t27 = (c5Eq*epsilon*t17*t18*t20*t21)/2.0;
  t28 = c3Eq+t26;
  t29 = c5Eq+t27;
  t0 = -coeff_hrxn*poro*(hRxnRate15+hRxnRate3*t29+hRxnRate1*(c4Eq+(c4Eq*epsilon*t17*t18*t20*t21)/2.0)-hRxnRate2*t28*t29-hRxnRate4*t28*(c2Eq+(c2Eq*epsilon*t17*t18*t20*t21)/2.0)-hRxnRate16*t28*(c7Eq+(c7Eq*epsilon*t17*t18*t20*t21)/2.0))+D1*coeff_diff*poro*t11*((c2Eq*epsilon*t2*t7*t13*t18*t19*t20)/2.0+c2Eq*epsilon*t2*t5*t18*t19*t20*t21*2.0+c2Eq*epsilon*t2*t9*t18*t19*t20*t21*2.0)-(D3*coeff_migr*e*poro*t10*t11*z3*((c3Eq*phiEq*t2*t3*t7*t22*t23*t24*pow(sin(t12),2.0))/4.0+c3Eq*phiEq*t2*t3*t9*t22*t24*t25*pow(sin(t15),2.0)-(epsilon*phiEq*t2*t7*t13*t17*t18*t20*t28)/2.0-epsilon*phiEq*t2*t5*t17*t18*t20*t21*t28*2.0-epsilon*phiEq*t2*t9*t17*t18*t20*t21*t28*2.0+c3Eq*phiEq*t2*t3*t5*(t19*t19)*t23*t24*t25))/kb+c3Eq*coeff_time*epsilon*poro*t10*t17*t18*t21*3.141592653589793*cos(t16);
