clear variables

% define basic symbolic variables 
syms x y z time

% reaction rates
syms hRxnRate1 hRxnRate2 hRxnRate3 hRxnRate4 hRxnRate5 hRxnRate6 hRxnRate7 hRxnRate8
syms hRxnRate9 hRxnRate10 hRxnRate11 hRxnRate12 hRxnRate13 hRxnRate14 hRxnRate15 hRxnRate16
syms fRxn1MoleReactant fRxn1MoleProduct1 fRxn1MoleProduct2
syms fRxn2MoleReactant fRxn2MoleProduct1 fRxn2MoleProduct2
syms fRxn1ExchangeCurrent fRxn1SurfaceRoughness fRxn1ActivationEnergy fRxn1TransferCoefficient fRxn1StandardElectrodePotential fRxn1CRef
syms fRxn2ExchangeCurrent fRxn2SurfaceRoughness fRxn2ActivationEnergy fRxn2TransferCoefficient fRxn2StandardElectrodePotential fRxn2CRef

% switch parameters, to turn effects on and off
syms coeff_time coeff_diff coeff_migr coeff_frxn coeff_hrxn

% define fields in term of other variables
syms c1(x,y,z,time)
syms c2(x,y,z,time)
syms c3(x,y,z,time)
syms c4(x,y,z,time)
syms c5(x,y,z,time)
syms c6(x,y,z,time)
syms c7(x,y,z,time)
syms phi(x,y,z,time)

% Equilibrium values 
syms c1Eq c2Eq c3Eq c4Eq c5Eq c6Eq c7Eq phiEq epsilon

% Domain size and problem parameters
syms Lx Ly Lz Period phiElectrode

% Physical constants
syms D1 D2 D3 D4 D5 D6 D7 poro tort na kb T e z1 z2 z3 z4 z5 z6 z7

% Manufactured solution fields
c1 = c1Eq + (epsilon * c1Eq * cos(2*pi*z/Lz) * sin(2*pi*x/Lx) * (1 + cos(pi*y/Ly)) / 2) * sin(2*pi*time/Period);
c2 = c2Eq + (epsilon * c2Eq * cos(2*pi*z/Lz) * cos(2*pi*x/Lx) * (1 + cos(pi*y/Ly)) / 2) * sin(2*pi*time/Period);
c3 = c3Eq + (epsilon * c3Eq * cos(2*pi*z/Lz) * cos(2*pi*x/Lx) * (1 + cos(pi*y/Ly)) / 2) * sin(2*pi*time/Period);
c4 = c4Eq + (epsilon * c4Eq * cos(2*pi*z/Lz) * cos(2*pi*x/Lx) * (1 + cos(pi*y/Ly)) / 2) * sin(2*pi*time/Period);
c5 = c5Eq + (epsilon * c5Eq * cos(2*pi*z/Lz) * cos(2*pi*x/Lx) * (1 + cos(pi*y/Ly)) / 2) * sin(2*pi*time/Period);
c6 = c6Eq + (epsilon * c6Eq * cos(2*pi*z/Lz) * cos(2*pi*x/Lx) * (1 + cos(pi*y/Ly)) / 2) * sin(2*pi*time/Period);
c7 = c7Eq + (epsilon * c7Eq * cos(2*pi*z/Lz) * cos(2*pi*x/Lx) * (1 + cos(pi*y/Ly)) / 2) * sin(2*pi*time/Period);
phi = phiEq + (epsilon * phiEq * cos(2*pi*z/Lz) * cos(2*pi*x/Lx) * (1 + cos(pi*y/Ly)) / 2) * sin(2*pi*time/Period);

% Temporal term
temp1 = poro * diff(c1, time);
temp2 = poro * diff(c2, time);
temp3 = poro * diff(c3, time);
temp4 = poro * diff(c4, time);
temp5 = poro * diff(c5, time);
temp6 = poro * diff(c6, time);
temp7 = poro * diff(c7, time);

% Diffusion term
diff1 = poro/tort * D1 * (diff(diff(c1, x), x) + diff(diff(c1, y), y) + diff(diff(c1, z), z));
diff2 = poro/tort * D2 * (diff(diff(c2, x), x) + diff(diff(c2, y), y) + diff(diff(c2, z), z));
diff3 = poro/tort * D3 * (diff(diff(c3, x), x) + diff(diff(c3, y), y) + diff(diff(c3, z), z));
diff4 = poro/tort * D4 * (diff(diff(c4, x), x) + diff(diff(c4, y), y) + diff(diff(c4, z), z));
diff5 = poro/tort * D5 * (diff(diff(c5, x), x) + diff(diff(c5, y), y) + diff(diff(c5, z), z));
diff6 = poro/tort * D6 * (diff(diff(c6, x), x) + diff(diff(c6, y), y) + diff(diff(c6, z), z));
diff7 = poro/tort * D7 * (diff(diff(c7, x), x) + diff(diff(c7, y), y) + diff(diff(c7, z), z));

% Electromigration term
migr1 = poro/tort * D1 * (z1*e/(kb*T)) * (diff(c1 * diff(phi, x), x) + diff(c1 * diff(phi, y), y) + diff(c1 * diff(phi, z), z));
migr2 = poro/tort * D2 * (z2*e/(kb*T)) * (diff(c2 * diff(phi, x), x) + diff(c2 * diff(phi, y), y) + diff(c2 * diff(phi, z), z));
migr3 = poro/tort * D3 * (z3*e/(kb*T)) * (diff(c3 * diff(phi, x), x) + diff(c3 * diff(phi, y), y) + diff(c3 * diff(phi, z), z));
migr4 = poro/tort * D4 * (z4*e/(kb*T)) * (diff(c4 * diff(phi, x), x) + diff(c4 * diff(phi, y), y) + diff(c4 * diff(phi, z), z));
migr5 = poro/tort * D5 * (z5*e/(kb*T)) * (diff(c5 * diff(phi, x), x) + diff(c5 * diff(phi, y), y) + diff(c5 * diff(phi, z), z));
migr6 = poro/tort * D6 * (z6*e/(kb*T)) * (diff(c6 * diff(phi, x), x) + diff(c6 * diff(phi, y), y) + diff(c6 * diff(phi, z), z));
migr7 = poro/tort * D7 * (z7*e/(kb*T)) * (diff(c7 * diff(phi, x), x) + diff(c7 * diff(phi, y), y) + diff(c7 * diff(phi, z), z));

% Homogeneous reactions
hRxn1 = poro * (-hRxnRate5 * c1 * c7 + hRxnRate6 * c5 - hRxnRate7 * c1 + hRxnRate8 * c4);
hRxn2 = poro * (hRxnRate3 * c5 - hRxnRate4 * c2 * c3 + hRxnRate9 * c5 * c5 - hRxnRate10 * c2 * c4 + hRxnRate13 * c5 * c7 - hRxnRate14 * c2);
hRxn3 = poro * (hRxnRate1 * c4 - hRxnRate2 * c3 * c5 + hRxnRate3 * c5 - hRxnRate4 * c2 * c3 + hRxnRate15 - hRxnRate16 * c3 * c7);
hRxn4 = poro * (-hRxnRate1 * c4 + hRxnRate2 * c5 * c3 + hRxnRate7 * c1 - hRxnRate8 * c4 + hRxnRate9 * c5 * c5 - hRxnRate10 * c2 * c4 - hRxnRate11 * c4 * c7 + hRxnRate12 * c5);
hRxn5 = poro * (hRxnRate1 * c4 - hRxnRate2 * c5 * c3 - hRxnRate3 * c5 + hRxnRate4 * c2 * c3 + hRxnRate5 * c1 * c7 - hRxnRate6 * c5 - hRxnRate9 * c5 * c5 + hRxnRate10 * c2 * c4 + hRxnRate11 * c4 * c7 - hRxnRate12 * c5 - hRxnRate13 * c5 * c7 + hRxnRate14 * c2);
hRxn6 = poro * (0);
hRxn7 = poro * (-hRxnRate5 * c1 * c7 + hRxnRate6 * c5 - hRxnRate11 * c4 * c7 + hRxnRate12 * c5 - hRxnRate13 * c5 * c7 + hRxnRate14 * c2 + hRxnRate15 - hRxnRate16 * c3 * c7);

% Faradaic reactions
fRxn1 = 1/(na * e) * fRxn1MoleReactant * fRxn1ExchangeCurrent * fRxn1SurfaceRoughness * exp(-fRxn1ActivationEnergy/(na * kb * T)) * exp(fRxn1TransferCoefficient * e/(kb*T) * (phiElectrode - phi - fRxn1StandardElectrodePotential)) * c1/fRxn1CRef;
fRxn7 = -1/(na * e) * fRxn1MoleProduct2 * fRxn1ExchangeCurrent * fRxn1SurfaceRoughness * exp(-fRxn1ActivationEnergy/(na * kb * T)) * exp(fRxn1TransferCoefficient * e/(kb*T) * (phiElectrode - phi - fRxn1StandardElectrodePotential)) * c1/fRxn1CRef ...
    - 1/(na * e) * fRxn2MoleProduct2 * fRxn2ExchangeCurrent * fRxn2SurfaceRoughness * exp(-fRxn2ActivationEnergy/(na * kb * T)) * exp(fRxn2TransferCoefficient * e/(kb*T) * (phiElectrode - phi - fRxn2StandardElectrodePotential));

% MMS Forcing Functions
force_c1 = simplify(coeff_time * temp1 - coeff_diff * diff1 - coeff_migr * migr1 - coeff_hrxn * hRxn1 - coeff_frxn * fRxn1);
force_c2 = simplify(coeff_time * temp2 - coeff_diff * diff2 - coeff_migr * migr2 - coeff_hrxn * hRxn2);
force_c3 = simplify(coeff_time * temp3 - coeff_diff * diff1 - coeff_migr * migr3 - coeff_hrxn * hRxn3);
force_c4 = simplify(coeff_time * temp4 - coeff_diff * diff1 - coeff_migr * migr4 - coeff_hrxn * hRxn4);
force_c5 = simplify(coeff_time * temp5 - coeff_diff * diff1 - coeff_migr * migr5 - coeff_hrxn * hRxn5);
force_c6 = simplify(coeff_time * temp6 - coeff_diff * diff6 - coeff_migr * migr6 - coeff_hrxn * hRxn6);
force_c7 = simplify(coeff_time * temp7 - coeff_diff * diff7 - coeff_migr * migr7 - coeff_hrxn * hRxn7 - coeff_frxn * fRxn7);
force_phi = simplify(z1 * e * (-coeff_diff * diff1 - coeff_migr * migr1 - coeff_frxn * fRxn1) ...
    + z2 * e * (-coeff_diff * diff2 - coeff_migr * migr2) ...
    + z3 * e * (-coeff_diff * diff3 - coeff_migr * migr3) ...
    + z4 * e * (-coeff_diff * diff4 - coeff_migr * migr4) ...
    + z5 * e * (-coeff_diff * diff5 - coeff_migr * migr5) ...
    + z6 * e * (-coeff_diff * diff6 - coeff_migr * migr6) ...
    + z7 * e * (-coeff_diff * diff7 - coeff_migr * migr7 - coeff_frxn * fRxn7));

args1 = symvar(force_c1);
args2 = symvar(force_c2);
args3 = symvar(force_c3);
args4 = symvar(force_c4);
args5 = symvar(force_c5);
args6 = symvar(force_c6);
args7 = symvar(force_c7);
argsphi = symvar(force_phi);

force_c1_f = matlabFunction(force_c1, 'Vars', args1);
force_c2_f = matlabFunction(force_c2, 'Vars', args2);
force_c3_f = matlabFunction(force_c3, 'Vars', args3);
force_c4_f = matlabFunction(force_c4, 'Vars', args4);
force_c5_f = matlabFunction(force_c5, 'Vars', args5);
force_c6_f = matlabFunction(force_c6, 'Vars', args6);
force_c7_f = matlabFunction(force_c7, 'Vars', args7);
force_phi_f = matlabFunction(force_phi, 'Vars', argsphi);

cargs1 = symvar(c1);
cargs2 = symvar(c2);
cargs3 = symvar(c3);
cargs4 = symvar(c4);
cargs5 = symvar(c5);
cargs6 = symvar(c6);
cargs7 = symvar(c7);
cargsphi = symvar(phi);

c1_f = matlabFunction(c1, 'Vars', cargs1);
c2_f = matlabFunction(c2, 'Vars', cargs2);
c3_f = matlabFunction(c3, 'Vars', cargs3);
c4_f = matlabFunction(c4, 'Vars', cargs4);
c5_f = matlabFunction(c5, 'Vars', cargs5);
c6_f = matlabFunction(c6, 'Vars', cargs6);
c7_f = matlabFunction(c7, 'Vars', cargs7);
phi_f = matlabFunction(phi, 'Vars', cargsphi);

save('mmsFunctions.mat', 'force_c1_f', 'force_c2_f', 'force_c3_f', 'force_c4_f', 'force_c5_f', 'force_c6_f', 'force_c7_f', 'force_phi_f', 'args1', 'args2', 'args3', 'args4', 'args5', 'args6', 'args7', 'argsphi', ...
    'c1_f', 'c2_f', 'c3_f', 'c4_f', 'c5_f', 'c6_f', 'c7_f', 'phi_f', 'cargs1', 'cargs2', 'cargs3', 'cargs4', 'cargs5', 'cargs6', 'cargs7', 'cargsphi')