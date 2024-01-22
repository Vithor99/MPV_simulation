% lambda_ss and Q_ss solver.%

function Q_ss=Q_solver(ups,sigma,eta, beta, iot, xi, varb, s, l_b_ss, l_g_ss, mu_e_ss, y_b, y_g, delta, sz,  EE_ss, UE_ss, EU_ss, AC_ss, U_ss)
options = optimoptions('fsolve','MaxIterations',100000, 'MaxFunctionEvaluations', 100000, 'FunctionTolerance',1.0000e-10, 'OptimalityTolerance', 1.0000e-10);
v0 = [0.5;0.5]; 
Fun = fsolve(@(x) [((beta/(1-beta*(1-delta)))*(iot*(((eta-1)/eta)*ups*(x(2)^((ups-1)/ups)))/(1+iot))*((1/(1+xi))*(((((eta-1)/eta)*ups*(x(2)^((ups-1)/ups)))*x(1)/varb)^xi))*mu_e_ss*((U_ss+delta*sz*(1-U_ss))*UE_ss*(((s*AC_ss*(1-delta)*UE_ss)/(EE_ss-delta+EU_ss))^iot)+(1-U_ss)*((EE_ss-delta*sz*UE_ss)/(AC_ss)))+(x(1)^(-sigma))-x(2));
(((1/(1+xi))*(((((eta-1)/eta)*ups*(x(2)^((ups-1)/ups)))*x(1)/varb)^xi))*(1+xi)*((y_g^(1+xi))*l_g_ss+(y_b^(1+xi))*l_b_ss)-(x(2)^(1/ups))-((x(2)-(x(1)^(-sigma)))/(iot*(((eta-1)/eta)*ups*(x(2)^((ups-1)/ups))))))], v0, options);
Q_ss = Fun(2);