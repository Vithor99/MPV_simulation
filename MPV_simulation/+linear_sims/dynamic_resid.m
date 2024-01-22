function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = linear_sims.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(63, 1);
residual(1) = y(83)-y(43)+y(44)-y(84);
residual(2) = y(46)-(1-params(17)*(1-params(16)))*((1+params(10))*(y(47)+y(57))+y(43)*params(10))-params(17)*(1-params(16))*(y(83)-y(43)+y(85));
residual(3) = T(1)*(y(43)+y(57)+y(47)+params(7)*y(48))-params(31)/params(46)*(y(83)+y(85));
residual(4) = y(43)+y(57)+y(47)+params(7)*y(49)-y(85)-y(83)-y(3);
lhs = y(50);
rhs = y(53)-y(55);
residual(5) = lhs - rhs;
residual(6) = T(3)-params(24)*params(47)*(params(48)+params(16)*(1-params(48))*params(22)+(1-params(16))*(1-params(48))*params(1))*y(51)-params(24)*params(47)*params(48)*(1-params(16)*params(22)-(1-params(16))*params(1))*y(4);
residual(7) = y(53)-(1-params(16))*(1-params(45)*params(1)*params(25)*params(4))*y(5)+y(49)*params(4)*params(45)*(1-params(16))*params(1)*params(25)-y(48)*params(44)*params(25)*(params(48)+params(16)*(1-params(48))*params(22))/params(49)*params(18)-y(4)*params(48)*params(18)*params(44)*(1-params(16)*params(22))*params(25)/params(49)+T(4)*(y(51)*params(9)+y(58));
residual(8) = params(50)*y(54)+y(53)*params(49)+params(48)*y(52);
residual(9) = y(52)-y(4)*(1-params(16)-params(44)*(1-params(16)*params(22))*params(25))+T(5)*(y(48)+y(51)*params(9)+y(58));
lhs = params(60)*(1+y(55));
rhs = 1-params(48)*(1+y(52));
residual(10) = lhs - rhs;
lhs = y(45);
rhs = y(84)*params(17)+T(6)*(y(47)+y(42)*(1-params(5))/params(5));
residual(11) = lhs - rhs;
residual(12) = y(51)*(1-params(9))-y(58)-y(47)-y(57)-(1+params(7))*(y(48)*T(8)+y(49)*(1-params(16))*(1-params(48))*params(1))/((1-params(16))*(1-params(48))*params(1)+T(8))+y(4)*T(9);
residual(13) = T(10)*(params(10)*(params(49)*T(11)+params(50)*T(12))*(y(57)+y(43)+y(47))-y(4)*params(48)*T(12)-y(5)*params(49)*(T(12)-T(11)))-T(13)*(y(42)*1/params(5)-y(57))-T(3)/(params(7)*params(43))-(y(47)+y(57))*(params(40)-T(2))/(params(7)*params(43));
residual(14) = y(44)-params(26)*y(1)-(1-params(26))*(params(27)*y(60)+y(52)*params(48)*params(28)+params(30)*params(56)*y(68))-y(59);
residual(15) = y(57)-params(14)*y(7)+x(it_, 3);
residual(16) = y(56)-params(12)*y(6)+x(it_, 1);
residual(17) = y(58)-params(15)*y(8)-x(it_, 4);
residual(18) = y(59)-params(13)*y(9)-x(it_, 2);
lhs = y(60);
rhs = y(45)+y(2)+y(10)+y(11)+y(12)+y(13)+y(14)+y(15)+y(16)+y(17)+y(18)+y(19);
residual(19) = lhs - rhs;
lhs = y(61);
rhs = y(51)*params(9)+y(58)+y(48)*T(14)+(1-T(14))*(y(49)+y(4)*params(48)/(1-params(48)))+params(45)*(1-params(16))*params(1)*params(25)/params(51)*params(4)*y(5)*params(49)/(1-params(48));
residual(20) = lhs - rhs;
lhs = y(62);
rhs = y(48)+y(51)*params(9)+y(58);
residual(21) = lhs - rhs;
lhs = y(63);
rhs = y(61)-y(62);
residual(22) = lhs - rhs;
lhs = y(64);
rhs = y(56)-y(43)*params(6);
residual(23) = lhs - rhs;
lhs = y(65);
rhs = y(51)+y(4)*T(15);
residual(24) = lhs - rhs;
lhs = y(66);
rhs = y(47)+y(57);
residual(25) = lhs - rhs;
lhs = y(67);
rhs = y(43)+y(46);
residual(26) = lhs - rhs;
lhs = y(68);
rhs = y(53)-y(55);
residual(27) = lhs - rhs;
lhs = y(69);
rhs = params(10)*(y(57)+y(43)+y(47));
residual(28) = lhs - rhs;
lhs = y(70);
rhs = y(42)-y(55)-y(69);
residual(29) = lhs - rhs;
lhs = y(71);
rhs = y(55)+y(61);
residual(30) = lhs - rhs;
lhs = y(72);
rhs = y(52)+y(62);
residual(31) = lhs - rhs;
lhs = y(20);
rhs = params(40)*exp(y(42));
residual(32) = lhs - rhs;
lhs = y(21);
rhs = params(41)*exp(y(43));
residual(33) = lhs - rhs;
lhs = y(22);
rhs = exp(y(44))*params(42);
residual(34) = lhs - rhs;
lhs = y(23);
rhs = params(43)*exp(y(47));
residual(35) = lhs - rhs;
lhs = y(34);
rhs = exp(y(57));
residual(36) = lhs - rhs;
lhs = y(24);
rhs = params(44)*exp(y(48));
residual(37) = lhs - rhs;
lhs = y(25);
rhs = params(45)*exp(y(49));
residual(38) = lhs - rhs;
lhs = y(31);
rhs = exp(y(56));
residual(39) = lhs - rhs;
lhs = y(26);
rhs = params(47)*exp(y(51));
residual(40) = lhs - rhs;
lhs = y(27);
rhs = params(48)*exp(y(52));
residual(41) = lhs - rhs;
lhs = y(28);
rhs = params(49)*exp(y(53));
residual(42) = lhs - rhs;
lhs = y(29);
rhs = params(50)*exp(y(54));
residual(43) = lhs - rhs;
lhs = y(32);
rhs = exp(y(58));
residual(44) = lhs - rhs;
lhs = y(33);
rhs = exp(y(59));
residual(45) = lhs - rhs;
lhs = y(35);
rhs = params(51)*exp(y(61));
residual(46) = lhs - rhs;
lhs = y(36);
rhs = exp(y(62))*params(52);
residual(47) = lhs - rhs;
lhs = y(37);
rhs = exp(y(63))*params(54);
residual(48) = lhs - rhs;
lhs = y(38);
rhs = exp(y(64))*params(58);
residual(49) = lhs - rhs;
lhs = y(39);
rhs = exp(y(65))*params(57);
residual(50) = lhs - rhs;
lhs = y(40);
rhs = exp(y(66))*params(59);
residual(51) = lhs - rhs;
lhs = y(30);
rhs = params(60)*exp(y(55));
residual(52) = lhs - rhs;
lhs = y(41);
rhs = params(17)*params(20)/((1-params(17)*(1-params(16)))*y(21));
residual(53) = lhs - rhs;
lhs = y(73);
rhs = y(2);
residual(54) = lhs - rhs;
lhs = y(74);
rhs = y(10);
residual(55) = lhs - rhs;
lhs = y(75);
rhs = y(11);
residual(56) = lhs - rhs;
lhs = y(76);
rhs = y(12);
residual(57) = lhs - rhs;
lhs = y(77);
rhs = y(13);
residual(58) = lhs - rhs;
lhs = y(78);
rhs = y(14);
residual(59) = lhs - rhs;
lhs = y(79);
rhs = y(15);
residual(60) = lhs - rhs;
lhs = y(80);
rhs = y(16);
residual(61) = lhs - rhs;
lhs = y(81);
rhs = y(17);
residual(62) = lhs - rhs;
lhs = y(82);
rhs = y(18);
residual(63) = lhs - rhs;

end
