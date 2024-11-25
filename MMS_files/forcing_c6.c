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
  t19 = sin(t16);
  t20 = t13+1.0;
  t21 = t17*t17;
  t22 = t18*t18;
  t23 = t19*t19;
  t24 = t20*t20;
  t25 = (c6Eq*epsilon*t17*t18*t19*t20)/2.0;
  t26 = c6Eq+t25;
  t0 = D6*coeff_diff*poro*t11*((c6Eq*epsilon*t2*t7*t13*t17*t18*t19)/2.0+c6Eq*epsilon*t2*t5*t17*t18*t19*t20*2.0+c6Eq*epsilon*t2*t9*t17*t18*t19*t20*2.0)-(D6*coeff_migr*e*poro*t10*t11*z6*((c6Eq*phiEq*t2*t3*t7*t21*t22*t23*pow(sin(t12),2.0))/4.0+c6Eq*phiEq*t2*t3*t5*t22*t23*t24*pow(sin(t14),2.0)+c6Eq*phiEq*t2*t3*t9*t21*t23*t24*pow(sin(t15),2.0)-(epsilon*phiEq*t2*t7*t13*t17*t18*t19*t26)/2.0-epsilon*phiEq*t2*t5*t17*t18*t19*t20*t26*2.0-epsilon*phiEq*t2*t9*t17*t18*t19*t20*t26*2.0))/kb+c6Eq*coeff_time*epsilon*poro*t10*t17*t18*t20*3.141592653589793*cos(t16);
