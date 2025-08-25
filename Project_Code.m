clear ALL
clc

% Constant Values
A = 37.5 * 10^(-4);
CMC = 50*10^(-3);
Rm = 21 * 10^(12);
t1 = 300;
Rfm = 7.6 * 10^(12);
K1 = 1.08 * 10^(-2);
Kc = 2 * 10^(5);
Kcv = 27;
S = 0.855;
Mu = 10^(-3);

% Guess Values
P = 413685;
Alpha = 10 * 10^(15);
Rv = 0.99;

% Case 1 - t<t1 .....................................

% Initial Values(At t = 0)
Vo = 0.1 * 10^(-3);
CoCV = 30 * 10^(-3);
CoRHL = 1;

syms CrCV(t3) CrRHL(t3) Vr(t3)

Cp = CrCV/(1+Kcv*S);
Rf = Rfm*(1-exp(-K1 * t3));

Vw = P/(Mu*(Rm+Rf));

ode1 = Vr * diff(CrCV,t3) == (Vw*A*Rv*(CrCV - Cp));
ode2 = diff(Vr,t3) == -Vw * A;
odes3 = [ode1; ode2];
S = dsolve(odes3);

Cond1 = Vr(0) == Vo;
Cond2 = CrCV(0) == CoCV;

Conds3 = [Cond1;Cond2];

Ssol = dsolve(odes3,Conds3);

Vrsol = Ssol.Vr;
CrCVsol = Ssol.CrCV;

ode3 = Vr * diff(CrRHL,t3) == (Vw*A*Rv*(CrRHL - CMC));

odes4 = [ode3;ode2];

S4 = dsolve(odes4);

Cond3 = CrRHL(0) == CoRHL;

Conds4 = [Cond1;Cond3];

Ssol3 = dsolve(odes4,Conds4);

CrRHLsol = Ssol3.CrRHL;

t3 = linspace(0,300,31);

Rf = Rfm*(1-exp(-K1 * t3));

CrCVvalue = subs(CrCVsol,{t3});
ValCrCV = double(CrCVvalue);
FinValCrCV = real(ValCrCV);

CrRHLvalue = subs(CrRHLsol,{t3});
ValCrRHL = double(CrRHLvalue);
FinValCrRHL = real(ValCrRHL);

VrValue = subs(Vrsol,{t3});
FinValVr = double(VrValue);

VwValue = subs(Vw,{t3});
VwVal1 = double(VwValue);

% Case 2 - t > t1.........................................

% Initial values (At t = 300)
Vo2 = FinValVr(end);
CoCV2 = FinValCrCV(end);
CoRHL2 = FinValCrRHL(end);
Vwo =  VwVal1(end);
Rfn = Rf(end);

syms CrCV(t4) CrRHL(t4) Vr(t4) Rc(t4)

S = 0.855;

Cp = CrCV/(1+Kcv*S);

Rc(t4) = (Alpha*(CoRHL2-CMC)*(Vo2 - Vr(t4)))/A;

VwN = P/(Mu*(Rm + Rfn + Rc(t4)));

ode4 = Vr * diff(CrCV,t4) == (VwN*A*Rv*(CrCV - Cp));
ode5 = diff(Vr,t4) == -VwN * A;
odes1 = [ode4;ode5];
S1 = dsolve(odes1);

Cond4 = Vr(300) == Vo2;
Cond5 = CrCV(300) == CoCV2;
Conds1 = [Cond4;Cond5];

Ssol2 = dsolve(odes1,Conds1);

Vrsol = Ssol2.Vr;
CrCVsol = Ssol2.CrCV;

ode6 = Vr * diff(CrRHL,t4) == (VwN*A*Rv*(CrRHL - CMC));
odes2 = [ode5;ode6];
S2 = dsolve(odes2);

Cond6 = Vr(300) == Vo2;
Cond7 = CrRHL(300) == CoRHL2;
Conds2 = [Cond6;Cond7];

Ssol3 = dsolve(odes2,Conds2);

Vrsol2 = Ssol3.Vr;
CrRHLsol = Ssol3.CrRHL;

t4 = linspace(300,3600,34);

CrCVvalue2 = subs(CrCVsol,{t4});
ValCrCV = double(CrCVvalue2);
FinValCrCV2 = real(ValCrCV);

CrRHLvalue2 = subs(CrRHLsol,{t4});
ValCrRHL = double(CrRHLvalue2);
FinValCrRHL2 = real(ValCrRHL);

VrValue2 = subs(Vrsol,{t4});
FinValVr2 = double(VrValue2);

Rc = (Alpha*(CoRHL2-CMC)*(Vo2 - FinValVr2))/A;

VwN = P./(Mu*(Rm + Rfn + Rc));

Rf = Rfm*(1-exp(-K1 * t3));

t3(end) = [];
FinValVr(end) = [];
FinValCrCV(end) = [];
FinValCrRHL(end) = [];
VwVal1(end)= [];

t = [t3 t4];

Vr = [FinValVr FinValVr2];
CrCV = [FinValCrCV FinValCrCV2];
CrRHL = [FinValCrRHL FinValCrRHL2];
Vw = [VwVal1 VwN];



