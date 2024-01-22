%% Moscarini Postel-Vinay log linearized around steady state
var
Q %MC final
lambda %Euler
R %taylor rule
X %MC service
r_u %intratemporal
r_e %intretemporal
teta %FEC
U %dynamic
l_b %dynamic
l_g %dynamic
empl %Employment 

mu %AR
phi_0 %AR
mps %AR
z %AR

EE % EE probability 
UE % UE probability 
AC
C
v
varphi
VU %Value of Unemployment (paychek of the unemployed)

%%%%%%%%%%%%%%%%%%%%%%% dev from ss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_hat %MC final
lambda_hat %Euler
R_hat %Taylor rule
pi %NKPC
W_hat %Presnt value
X_hat %MC service
r_u_hat %intratemporal
r_e_hat %intretemporal
mu_e_hat 
teta_hat %FEC
U_hat %dynamic
l_b_hat %dynamic
l_g_hat %dynamic
empl_hat

mu_hat %AR
z_hat %AR
phi_0_hat %AR
mps_hat %AR


INFL
EE_hat
UE_hat
AC_hat
C_hat
v_hat
varphi_hat
WplusLam %i.e. W adjusted by utility  
b_share %share of bad matches from one period before 
H_hat %total hours of work
TFP_hat %Total Factor productivity Q/(empl*H)
EEf_hat %EE flow
UEf_hat %UE flow
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo eps_mu %innovation on preference for consumption
eps_mps %innovation on monetary policy shock
eps_z %innovation on tech shock 
eps_phi0 %innovation on matching technology
;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameters
%we start with our varying parameters
s      %param1 %search intensity of workers
y_b    %param 2 %productivity of good match
s_g    %efficiency gap between bad and good matches
xi_g   %param 3 % cond. on matching, prob that match is good type 
ups    %Returns to scale 
sigma  %intert elasticity of substitution
iot    %recruiting cost elasticity                                                     
varb   %weight on cost of effort (intensive margin)
psi    %meeting function elasticity                                                        
xi     %Frisch elasticity of effort                                                      
zeta   %Calvo Parameter

%all the rest
rho_mu % AR coefficient on preference shock
rho_mps % AR coefficient on monetary policy shock
rho_z % AR coefficient on technology shock
rho_phi0 % AR coefficient on matching technology shock

delta % probability of job destruction
beta %discount factor
xi_b % cond. on matching, prob that match is bad type (1-xi_g)
y_g %productivity of bad match 
b %flow value of leisure (extensive margin)
eta %elasticity of substitution for consumption
sz  %unemployed search efficiency at short duration
kap_s %recruting cost
c_p %cost of vacancy
phi_cap_ss
%for Taylor rule
rho_r %monetary rule coeff.1 (degree of IR smoothing)
phi_pi %monetary rule coeff.2 
phi_u %monetary rule coeff.3
phi_EE
phi_b_share
mu_u
sigma_z
sigma_mps
sigma_mu
sigma_phi0

b_cap
varb_cap
prox_u
prox_e

%steady state values 
Q_ss %MC final
lambda_ss %Euler
R_ss %taylor rule
X_ss %MC service
r_u_ss %intratemporal
r_e_ss %intretemporal
mu_e_ss 
teta_ss %FEC
U_ss %dynamic
l_b_ss %dynamic
l_g_ss %dynamic

EE_ss
UE_ss
EU_ss
AC_ss
EEf_ss
b_share_ss
v_ss
C_ss
varphi_ss
empl_ss
;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters estimated via SMM in mod1
%initial values
s = 0.99;    %param 1  (0.5)                                                            
y_b= 1;  %param 2 
sg = 0.01; %0.8
xi_g= 0.5;  %param 3
ups = 0.67; 
sigma = 0.33;
iot=2.35;                                                              
varb=0.37;
psi=0.38; %alpha                                                             
xi= 0.37;                                                             
zeta = 0.11; 
 
% Empiricals (debatable) 
EE_ss= 0.024;
UE_ss=0.25;
EU_ss=0.01314;
AC_ss=EE_ss/UE_ss; %by definition
U_ss=EU_ss/(UE_ss+EU_ss); %UR consistent with EU and UE 
empl_ss=1-U_ss;
EEf_ss = EE_ss*empl_ss;

%By assumption set:
teta_ss=1; 
phi_cap_ss=1; %follows

%step 1: Precalibrated following Moscarini and Postel vinaj table 
eta=6;
beta = 0.995; %discount factor(to get r=0.005)
rho_mu = 0.94; 
rho_z = 0.494; 
rho_phi0=0.5; 
R_ss = (1/beta)-1; 

%shock variance
sigma_z=0.0367^2;
sigma_mps=0.023^2; 
sigma_mu=0.023^2;
sigma_phi0=0.023^2; %random

%step 2: estimate by GMM Taylor rules parameters
rho_mps = 0.93; % AR coefficient on monetary policy shock (Moscarini)
phi_pi = 1.17;  %monetary rule coeff.2 
phi_u = -0.05;  %monetary rule coeff.3
phi_EE = 0.05; %CHANGE
phi_b_share = 0.05; %CHANGE
rho_r = 0.85;   %monetary rule coeff.1 (degree of IR smoothing) %CHANGE

%step 3: Give parametrization for mathces quality and distribution 
y_g = y_b*(1+sg);
xi_b=1-xi_g;

%step 4: Compute delta and sz
delta = (xi_g*(EE_ss+EU_ss))/(xi_g*xi_b-AC_ss+xi_g);
sz=(1-(EU_ss/delta))/UE_ss; 


%step 5: r_ratio 
r_ratio = (s*AC_ss*(1-delta)*UE_ss)/(EE_ss-delta+EU_ss); 

%step 6: l_g and l_b
l_b_ss = (delta*(1-U_ss)*xi_b)/(delta+xi_g*((EE_ss-delta*sz*UE_ss)/(AC_ss)));
l_g_ss = 1- U_ss - l_b_ss; 
b_share_ss = l_b_ss/l_g_ss;

%step 7: mu_u and mu_e
mu_u=(y_g^(1+xi))*xi_g+(y_b^(1+xi))*xi_b;
mu_e_ss = (xi_g*(l_b_ss/(1-U_ss))*((y_g^(1+xi))-(y_b^(1+xi))));

%step 8: solve for X, lambda and Q (FUNCTION in steady.m file)
lambda_ss = lambda_solver(ups,sigma,eta, beta, iot, xi, varb, s, l_b_ss, l_g_ss, mu_e_ss, y_b, y_g, delta, sz,  EE_ss, UE_ss, EU_ss, AC_ss, U_ss);
Q_ss = Q_solver(ups,sigma,eta, beta, iot, xi, varb, s, l_b_ss, l_g_ss, mu_e_ss, y_b, y_g, delta, sz,  EE_ss, UE_ss, EU_ss, AC_ss, U_ss); 

X_ss = (((eta-1)/eta)*ups*(Q_ss^((ups-1)/ups)));
varb_cap = ((1/(1+xi))*((X_ss*lambda_ss/varb)^xi));
b_cap = varb_cap*mu_u-varb_cap*mu_e_ss*(((s*AC_ss*(1-delta)*UE_ss)/(EE_ss-delta+EU_ss))^iot);

%step 9: b
b = b_cap*X_ss*lambda_ss; 

%step 10: let prox1=kap_s*r_u_ss^iot and prox2=kap_s*r_e_ss^iot
prox_u = (beta/(1-beta*(1-delta)))*(varb_cap*mu_u-b_cap);
prox_e = (beta/(1-beta*(1-delta)))*varb_cap*mu_e_ss;

%step 11: c_p*teta which is equivalent to c_p since teta=1
c_p = (Q_ss-(lambda_ss^(-sigma)))/(U_ss+(delta*sz+(1-delta)*s)*(1-U_ss));

%final assumption
r_u_ss=UE_ss;
r_e_ss = (EE_ss-delta*sz*UE_ss)/(AC_ss*(1-delta)*s);
kap_s = prox_u / (r_u_ss^iot); 

%remaining variables
C_ss=lambda_ss^(-sigma);
v_ss= U_ss+(delta*sz+(1-delta)*s)*(1-U_ss);
varphi_ss = X_ss; 
VU_ss = (beta*b)/(lambda_ss*(1-beta*(1-delta)));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model; 

lambda_hat(+1)-lambda_hat+R_hat-pi(+1)=0; %Euler

W_hat - (1-beta*(1-delta))*((1+xi)*(X_hat+z_hat)+xi*lambda_hat)-beta*(1-delta)*(lambda_hat(+1)-lambda_hat+W_hat(+1))=0; %Present Value of Service relative price

((r_u_ss/r_e_ss)^(iot))*(iot*r_u_hat+X_hat+z_hat+lambda_hat)-(mu_u/mu_e_ss)*(W_hat(+1)+lambda_hat(+1))=0; %Recruting intensity of unemployed job applicant 

iot*r_e_hat+X_hat+z_hat-W_hat(+1)-lambda_hat(+1)+lambda_hat-mu_e_hat(-1)=0; %Recruting intensity of employed job applicant 

mu_e_hat = l_b_hat - empl_hat; %Surplus of good matches

Q_ss*Q_hat+(lambda_ss^(-sigma))*(sigma*lambda_hat-mu_hat)-c_p*teta_ss*(U_ss+delta*(1-U_ss)*sz+(1-delta)*(1-U_ss)*s)*teta_hat-c_p*teta_ss*U_ss*(1-delta*sz-(1-delta)*s)*U_hat(-1)=0; %final good market clearing 

l_b_hat-(1-delta)*(1-s*phi_cap_ss*r_e_ss*xi_g)*l_b_hat(-1)+(1-delta)*s*phi_cap_ss*r_e_ss*xi_g*r_e_hat-((U_ss+delta*(1-U_ss)*sz)/l_b_ss)*phi_cap_ss*r_u_ss*xi_b*r_u_hat-(((1-delta*sz)*phi_cap_ss*r_u_ss*xi_b)/(l_b_ss))*U_ss*U_hat(-1)+((1-delta)*s*r_e_ss*xi_g-((U_ss+delta*(1-U_ss)*sz)/(l_b_ss))*r_u_ss*xi_b)*phi_cap_ss*(psi*teta_hat+phi_0_hat)=0; %Bad matches dynamic

l_g_ss*l_g_hat+l_b_ss*l_b_hat+U_ss*U_hat=0; %good matches dynamic 

U_hat-(1-delta-(1-delta*sz)*phi_cap_ss*r_u_ss)*U_hat(-1)+(1-delta*sz+(delta*sz/U_ss))*phi_cap_ss*r_u_ss*(psi*teta_hat+phi_0_hat+r_u_hat)=0; %Unemployment dynamic

empl_ss*(1+empl_hat)=1-U_ss*(1+U_hat); 

#kappa = (zeta/(1-zeta))*((1-beta*(1-zeta))/(1+(eta/ups)-eta));
pi = kappa*(((1-ups)/ups)*Q_hat + X_hat) + beta * pi(+1);               %NKPC in logs

(1-psi)*teta_hat-phi_0_hat-X_hat-z_hat-(1+iot)*(((U_ss+delta*(1-U_ss)*sz)*((r_u_ss/r_e_ss)^(1+iot))*r_u_hat + (1-delta)*(1-U_ss)*s*r_e_hat)/((U_ss+delta*(1-U_ss)*sz)*((r_u_ss/r_e_ss)^(1+iot))+(1-delta)*(1-U_ss)*s))+(((1-delta*sz-(1-delta)*s)/(U_ss+(delta*sz+(1-delta)*s)*(1-U_ss)))-(((1-delta*sz)*((r_u_ss/r_e_ss)^(1+iot))-(1-delta)*s)/((U_ss+delta*(1-U_ss)*sz)*((r_u_ss/r_e_ss)^(1+iot))+(1-delta)*(1-U_ss)*s)))*U_ss*U_hat(-1)=0; 

((X_ss*lambda_ss/varb)^xi)*(xi*(l_b_ss*(y_b^(1+xi))+l_g_ss*(y_g^(1+xi)))*(lambda_hat+X_hat+z_hat)-(y_g^(1+xi))*U_ss*U_hat(-1)-((y_g^(1+xi))-(y_b^(1+xi)))*l_b_ss*l_b_hat(-1))-(Q_ss^(1/ups))*((1/ups)*Q_hat-z_hat)-((Q_ss*Q_hat+(lambda_ss^(-sigma))*(sigma*lambda_hat-mu_hat))/(iot*X_ss))-((Q_ss-(lambda_ss^(-sigma)))/(iot*X_ss))*(z_hat+X_hat)=0;

@#ifdef BenchmarkModel
    R_hat-rho_r*R_hat(-1)-(1-rho_r)*(phi_pi*INFL+phi_u*U_ss*U_hat)-mps_hat=0; %Taylor rule
@#endif

@#ifdef EE_SimpleRule
    R_hat-rho_r*R_hat(-1)-(1-rho_r)*(phi_pi*INFL+phi_u*U_ss*U_hat+ phi_EE * EEf_ss * EEf_hat)-mps_hat=0; %Taylor rule
@#endif

@#ifdef b_shareRule
    R_hat-rho_r*R_hat(-1)-(1-rho_r)*(phi_pi*INFL+phi_u*U_ss*U_hat+ phi_b_share * b_share_ss * b_share)-mps_hat=0; %Taylor rule
@#endif


z_hat-rho_z*z_hat(-1)+eps_z=0; %technology AR
mu_hat-rho_mu*mu_hat(-1)+eps_mu=0; %technology AR
phi_0_hat-rho_phi0*phi_0_hat(-1)-eps_phi0=0; %technology AR
mps_hat-rho_mps*mps_hat(-1)-eps_mps=0; %technology AR

INFL =pi+pi(-1)+pi(-2)+pi(-3)+pi(-4)+pi(-5)+pi(-6)+pi(-7)+pi(-8)+pi(-9)+pi(-10)+pi(-11);
EE_hat = phi_0_hat+psi*teta_hat+(delta*sz*phi_cap_ss*r_u_ss/EE_ss)*r_u_hat+(1-(delta*sz*phi_cap_ss*r_u_ss/EE_ss))*(r_e_hat+(U_ss/(1-U_ss))*U_hat(-1))+((1-delta)*s*phi_cap_ss*r_e_ss/EE_ss)*(xi_g*(l_b_ss*l_b_hat(-1)/(1-U_ss)));
UE_hat = psi*teta_hat+phi_0_hat+r_u_hat;
AC_hat=EE_hat-UE_hat; 
C_hat = mu_hat - sigma*lambda_hat; 
v_hat = teta_hat+((1-(1-delta)*s-delta*sz)/(U_ss+(delta*sz+(1-delta)*s)*(1-U_ss)))*U_ss*U_hat(-1);
varphi_hat=X_hat+z_hat; 
WplusLam = W_hat + lambda_hat; 
b_share = l_b_hat-empl_hat;
H_hat = xi*(lambda_hat+X_hat+z_hat);
TFP_hat = Q_hat - empl_hat - H_hat;
EEf_hat = empl_hat + EE_hat; 
UEf_hat = U_hat + UE_hat; 

%% non linear variables
Q = exp(Q_hat)*Q_ss;%MC final
lambda = exp(lambda_hat)*lambda_ss;%Euler
R=exp(R_hat)*R_ss; %taylor rule
X=exp(X_hat)*X_ss; %MC service
z=exp(z_hat); %AR
r_u=exp(r_u_hat)*r_u_ss; %intratemporal
r_e=exp(r_e_hat)*r_e_ss; %intretemporal 
mu=exp(mu_hat); %AR
teta=exp(teta_hat)*teta_ss; %FEC
U =exp(U_hat)*U_ss;%dynamic
l_b=exp(l_b_hat)*l_b_ss; %dynamic
l_g=exp(l_g_hat)*l_g_ss; %dynamic
phi_0=exp(phi_0_hat); %AR
mps=exp(mps_hat); %AR
EE=exp(EE_hat)*EE_ss; 
UE=exp(UE_hat)*UE_ss;
AC=exp(AC_hat)*AC_ss;
C=exp(C_hat)*C_ss; 
v=exp(v_hat)*v_ss;
varphi=exp(varphi_hat)*varphi_ss; 
empl=exp(empl_hat)*empl_ss; 
VU = beta*b/(lambda*(1-beta*(1-delta)));
end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steady_state_model;
Q=Q_ss; %MC final
lambda=lambda_ss; %Euler
R=R_ss; %taylor rule
X=X_ss; %MC service
z=1; %AR
r_u=r_u_ss; %intratemporal
r_e=r_e_ss; %intretemporal
mu =1;%AR
teta=teta_ss; %FEC
U=U_ss; %dynamic
l_b=l_b_ss; %dynamic
l_g=l_g_ss; %dynamic
phi_0=1; %AR
mps=1; %AR

%%%%%%%%%%%%%%%%%%%%%%% Auxiliary %%%%%%%%%%%%%%%%%%
EE=EE_ss;
UE=UE_ss;
AC=AC_ss;
C=C_ss;
v=v_ss;
varphi=varphi_ss; 
empl=empl_ss;
VU=beta*b/(lambda_ss*(1-beta*(1-delta)));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@#ifdef RamseyModel
    planner_objective INFL^2 +1*Q_hat^2;
    ramsey_model(instruments=(R_hat),planner_discount=beta);
@#endif

steady;
resid;
check;

shocks;    
var eps_mu = sigma_mu;
var eps_mps = sigma_mps; %(Moscarini) 
var eps_z = sigma_z; %Moscarini
end;


stoch_simul(order = 1
    ,irf=0
    ,periods=10000
    ,irf_plot_threshold=0) mps EEf_hat EE_hat U_hat C_hat b_share;







