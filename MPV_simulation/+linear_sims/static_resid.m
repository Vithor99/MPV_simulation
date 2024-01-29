function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = linear_sims.static_resid_tt(T, y, x, params);
end
residual = zeros(64, 1);
residual(1) = y(24)-y(25);
residual(2) = y(26)-(1-params(17)*(1-params(16)))*((1+params(10))*(y(27)+y(39))+y(23)*params(10))-y(26)*params(17)*(1-params(16));
residual(3) = T(2)*(y(23)+y(39)+y(27)+params(7)*y(28))-params(29)/params(44)*(y(23)+y(26));
residual(4) = y(39)+y(27)+params(7)*y(29)-y(26)-y(30);
lhs = y(30);
rhs = y(34)-y(36);
residual(5) = lhs - rhs;
residual(6) = T(4)-params(23)*params(45)*(params(46)+params(16)*(1-params(46))*params(21)+(1-params(16))*(1-params(46))*params(1))*y(31)-params(23)*params(45)*params(46)*(1-params(16)*params(21)-(1-params(16))*params(1))*y(33);
residual(7) = y(34)-y(34)*(1-params(16))*(1-params(43)*params(1)*params(24)*params(4))+y(29)*params(4)*params(43)*(1-params(16))*params(1)*params(24)-y(28)*params(42)*params(24)*(params(46)+params(16)*(1-params(46))*params(21))/params(47)*params(18)-y(33)*params(46)*params(18)*params(42)*(1-params(16)*params(21))*params(24)/params(47)+T(5)*y(31)*params(9);
residual(8) = params(48)*y(35)+y(34)*params(47)+params(46)*y(33);
residual(9) = y(33)-y(33)*(1-params(16)-params(42)*(1-params(16)*params(21))*params(24))+T(6)*(y(28)+y(31)*params(9));
lhs = params(59)*(1+y(36));
rhs = 1-params(46)*(1+y(33));
residual(10) = lhs - rhs;
lhs = y(25);
rhs = params(17)*y(25)+T(1)*(y(27)+y(22)*(1-params(5))/params(5));
residual(11) = lhs - rhs;
lhs = y(32);
rhs = y(31)*params(9);
residual(12) = lhs - rhs;
residual(13) = y(31)*(1-params(9))-y(27)-y(39)-(1+params(7))*(y(28)*T(8)+y(29)*(1-params(16))*(1-params(46))*params(1))/((1-params(16))*(1-params(46))*params(1)+T(8))+y(33)*T(9);
residual(14) = T(10)*(params(10)*(params(47)*T(11)+params(48)*T(12))*(y(39)+y(23)+y(27))-y(33)*params(46)*T(12)-y(34)*params(47)*(T(12)-T(11)))-T(13)*(y(22)*1/params(5)-y(39))-T(4)/(params(7)*params(41))-(y(27)+y(39))*(params(38)-T(3))/(params(7)*params(41));
residual(15) = y(24)-y(24)*params(25)-(1-params(25))*(y(25)*params(26)+y(22)*params(38)*params(27)+params(28)*params(52)*y(44))-y(40);
residual(16) = y(39)-y(39)*params(14)+x(2);
residual(17) = y(38)-y(38)*params(12)+x(1);
residual(18) = y(40)-y(40)*params(13)-x(3);
lhs = y(41);
rhs = y(25)+y(25)+y(25)+y(25)+y(25)+y(25)+y(25)+y(25)+y(25)+y(25)+y(25)+y(25);
residual(19) = lhs - rhs;
lhs = y(42);
rhs = y(31)*params(9)+y(28)*T(14)+(1-T(14))*(y(29)+y(33)*params(46)/(1-params(46)))+params(43)*(1-params(16))*params(1)*params(24)/params(49)*params(4)*y(34)*params(47)/(1-params(46));
residual(20) = lhs - rhs;
lhs = y(43);
rhs = y(28)+y(31)*params(9);
residual(21) = lhs - rhs;
lhs = y(44);
rhs = y(42)-y(43);
residual(22) = lhs - rhs;
lhs = y(45);
rhs = y(38)-y(23)*params(6);
residual(23) = lhs - rhs;
lhs = y(46);
rhs = y(31)+y(33)*T(15);
residual(24) = lhs - rhs;
lhs = y(47);
rhs = y(27)+y(39);
residual(25) = lhs - rhs;
lhs = y(48);
rhs = y(23)+y(26);
residual(26) = lhs - rhs;
lhs = y(49);
rhs = y(34)-y(36);
residual(27) = lhs - rhs;
lhs = y(50);
rhs = y(35)-y(36);
residual(28) = lhs - rhs;
lhs = y(51);
rhs = params(10)*(y(39)+y(23)+y(27));
residual(29) = lhs - rhs;
lhs = y(37);
rhs = y(36)+y(51);
residual(30) = lhs - rhs;
lhs = y(52);
rhs = y(22)-y(36)-y(51);
residual(31) = lhs - rhs;
lhs = y(53);
rhs = y(36)+y(42);
residual(32) = lhs - rhs;
lhs = y(54);
rhs = y(33)+y(43);
residual(33) = lhs - rhs;
lhs = y(1);
rhs = params(38)*exp(y(22));
residual(34) = lhs - rhs;
lhs = y(2);
rhs = params(39)*exp(y(23));
residual(35) = lhs - rhs;
lhs = y(3);
rhs = exp(y(24))*params(40);
residual(36) = lhs - rhs;
lhs = y(4);
rhs = params(41)*exp(y(27));
residual(37) = lhs - rhs;
lhs = y(14);
rhs = exp(y(39));
residual(38) = lhs - rhs;
lhs = y(5);
rhs = params(42)*exp(y(28));
residual(39) = lhs - rhs;
lhs = y(6);
rhs = params(43)*exp(y(29));
residual(40) = lhs - rhs;
lhs = y(12);
rhs = exp(y(38));
residual(41) = lhs - rhs;
lhs = y(7);
rhs = params(45)*exp(y(31));
residual(42) = lhs - rhs;
lhs = y(8);
rhs = params(46)*exp(y(33));
residual(43) = lhs - rhs;
lhs = y(9);
rhs = params(47)*exp(y(34));
residual(44) = lhs - rhs;
lhs = y(10);
rhs = params(48)*exp(y(35));
residual(45) = lhs - rhs;
lhs = y(13);
rhs = exp(y(40));
residual(46) = lhs - rhs;
lhs = y(15);
rhs = params(49)*exp(y(42));
residual(47) = lhs - rhs;
lhs = y(16);
rhs = exp(y(43))*params(50);
residual(48) = lhs - rhs;
lhs = y(17);
rhs = params(52)*exp(y(44));
residual(49) = lhs - rhs;
lhs = y(18);
rhs = exp(y(45))*params(57);
residual(50) = lhs - rhs;
lhs = y(19);
rhs = exp(y(46))*params(56);
residual(51) = lhs - rhs;
lhs = y(20);
rhs = exp(y(47))*params(58);
residual(52) = lhs - rhs;
lhs = y(11);
rhs = params(59)*exp(y(36));
residual(53) = lhs - rhs;
lhs = y(21);
rhs = params(17)*params(19)/((1-params(17)*(1-params(16)))*y(2));
residual(54) = lhs - rhs;
lhs = y(55);
rhs = y(25);
residual(55) = lhs - rhs;
lhs = y(56);
rhs = y(25);
residual(56) = lhs - rhs;
lhs = y(57);
rhs = y(25);
residual(57) = lhs - rhs;
lhs = y(58);
rhs = y(25);
residual(58) = lhs - rhs;
lhs = y(59);
rhs = y(25);
residual(59) = lhs - rhs;
lhs = y(60);
rhs = y(25);
residual(60) = lhs - rhs;
lhs = y(61);
rhs = y(25);
residual(61) = lhs - rhs;
lhs = y(62);
rhs = y(25);
residual(62) = lhs - rhs;
lhs = y(63);
rhs = y(25);
residual(63) = lhs - rhs;
lhs = y(64);
rhs = y(25);
residual(64) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
