function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = linear_sims.static_g1_tt(T, y, x, params);
end
g1 = zeros(63, 63);
g1(1,25)=1;
g1(1,26)=(-1);
g1(2,24)=(-((1-params(17)*(1-params(16)))*params(10)));
g1(2,27)=1-params(17)*(1-params(16));
g1(2,28)=(-((1-params(17)*(1-params(16)))*(1+params(10))));
g1(2,38)=(-((1-params(17)*(1-params(16)))*(1+params(10))));
g1(3,24)=T(2)-params(31)/params(46);
g1(3,27)=(-(params(31)/params(46)));
g1(3,28)=T(2);
g1(3,29)=params(7)*T(2);
g1(3,38)=T(2);
g1(4,27)=(-1);
g1(4,28)=1;
g1(4,30)=params(7);
g1(4,31)=(-1);
g1(4,38)=1;
g1(5,31)=1;
g1(5,34)=(-1);
g1(5,36)=1;
g1(6,23)=params(40);
g1(6,24)=params(6)*T(3);
g1(6,32)=(-(params(24)*params(47)*(params(48)+params(16)*(1-params(48))*params(22)+(1-params(16))*(1-params(48))*params(1))));
g1(6,33)=(-(params(24)*params(47)*params(48)*(1-params(16)*params(22)-(1-params(16))*params(1))));
g1(6,37)=(-T(3));
g1(7,29)=(-(params(44)*params(25)*(params(48)+params(16)*(1-params(48))*params(22))/params(49)*params(18)));
g1(7,30)=params(4)*params(45)*(1-params(16))*params(1)*params(25);
g1(7,32)=T(5)*params(9);
g1(7,33)=(-(params(48)*params(18)*params(44)*(1-params(16)*params(22))*params(25)/params(49)));
g1(7,34)=1-(1-params(16))*(1-params(45)*params(1)*params(25)*params(4));
g1(7,39)=T(5);
g1(8,33)=params(48);
g1(8,34)=params(49);
g1(8,35)=params(50);
g1(9,29)=T(6);
g1(9,32)=params(9)*T(6);
g1(9,33)=1-(1-params(16)-params(44)*(1-params(16)*params(22))*params(25));
g1(9,39)=T(6);
g1(10,33)=params(48);
g1(10,36)=params(60);
g1(11,23)=(-(T(1)*(1-params(5))/params(5)));
g1(11,26)=1-params(17);
g1(11,28)=(-T(1));
g1(12,28)=(-1);
g1(12,29)=(-((1+params(7))*T(8)/((1-params(16))*(1-params(48))*params(1)+T(8))));
g1(12,30)=(-((1+params(7))*(1-params(16))*(1-params(48))*params(1)/((1-params(16))*(1-params(48))*params(1)+T(8))));
g1(12,32)=1-params(9);
g1(12,33)=T(9);
g1(12,38)=(-1);
g1(12,39)=(-1);
g1(13,23)=(-(1/params(5)*T(13)))-params(40)/(params(7)*params(43));
g1(13,24)=T(10)*params(10)*(params(49)*T(11)+params(50)*T(12))-params(6)*T(3)/(params(7)*params(43));
g1(13,28)=T(10)*params(10)*(params(49)*T(11)+params(50)*T(12))-(params(40)-T(3))/(params(7)*params(43));
g1(13,33)=T(10)*(-(params(48)*T(12)));
g1(13,34)=T(10)*(-(params(49)*(T(12)-T(11))));
g1(13,37)=(-((-T(3))/(params(7)*params(43))));
g1(13,38)=T(13)+T(10)*params(10)*(params(49)*T(11)+params(50)*T(12))-(params(40)-T(3))/(params(7)*params(43));
g1(14,25)=1-params(26);
g1(14,33)=(-((1-params(26))*params(48)*params(28)));
g1(14,40)=(-1);
g1(14,41)=(-((1-params(26))*params(27)));
g1(14,49)=(-((1-params(26))*params(30)*params(56)));
g1(15,38)=1-params(14);
g1(16,37)=1-params(12);
g1(17,39)=1-params(15);
g1(18,40)=1-params(13);
g1(19,26)=(-12);
g1(19,41)=1;
g1(20,29)=(-T(14));
g1(20,30)=(-(1-T(14)));
g1(20,32)=(-params(9));
g1(20,33)=(-((1-T(14))*params(48)/(1-params(48))));
g1(20,34)=(-(params(45)*(1-params(16))*params(1)*params(25)/params(51)*params(4)*params(49)/(1-params(48))));
g1(20,39)=(-1);
g1(20,42)=1;
g1(21,29)=(-1);
g1(21,32)=(-params(9));
g1(21,39)=(-1);
g1(21,43)=1;
g1(22,42)=(-1);
g1(22,43)=1;
g1(22,44)=1;
g1(23,24)=params(6);
g1(23,37)=(-1);
g1(23,45)=1;
g1(24,32)=(-1);
g1(24,33)=(-T(15));
g1(24,46)=1;
g1(25,28)=(-1);
g1(25,38)=(-1);
g1(25,47)=1;
g1(26,24)=(-1);
g1(26,27)=(-1);
g1(26,48)=1;
g1(27,34)=(-1);
g1(27,36)=1;
g1(27,49)=1;
g1(28,24)=(-params(10));
g1(28,28)=(-params(10));
g1(28,38)=(-params(10));
g1(28,50)=1;
g1(29,23)=(-1);
g1(29,36)=1;
g1(29,50)=1;
g1(29,51)=1;
g1(30,36)=(-1);
g1(30,42)=(-1);
g1(30,52)=1;
g1(31,33)=(-1);
g1(31,43)=(-1);
g1(31,53)=1;
g1(32,1)=1;
g1(32,23)=(-(params(40)*exp(y(23))));
g1(33,2)=1;
g1(33,24)=(-(params(41)*exp(y(24))));
g1(34,3)=1;
g1(34,25)=(-(exp(y(25))*params(42)));
g1(35,4)=1;
g1(35,28)=(-(params(43)*exp(y(28))));
g1(36,15)=1;
g1(36,38)=(-exp(y(38)));
g1(37,5)=1;
g1(37,29)=(-(params(44)*exp(y(29))));
g1(38,6)=1;
g1(38,30)=(-(params(45)*exp(y(30))));
g1(39,12)=1;
g1(39,37)=(-exp(y(37)));
g1(40,7)=1;
g1(40,32)=(-(params(47)*exp(y(32))));
g1(41,8)=1;
g1(41,33)=(-(params(48)*exp(y(33))));
g1(42,9)=1;
g1(42,34)=(-(params(49)*exp(y(34))));
g1(43,10)=1;
g1(43,35)=(-(params(50)*exp(y(35))));
g1(44,13)=1;
g1(44,39)=(-exp(y(39)));
g1(45,14)=1;
g1(45,40)=(-exp(y(40)));
g1(46,16)=1;
g1(46,42)=(-(params(51)*exp(y(42))));
g1(47,17)=1;
g1(47,43)=(-(exp(y(43))*params(52)));
g1(48,18)=1;
g1(48,44)=(-(exp(y(44))*params(54)));
g1(49,19)=1;
g1(49,45)=(-(exp(y(45))*params(58)));
g1(50,20)=1;
g1(50,46)=(-(exp(y(46))*params(57)));
g1(51,21)=1;
g1(51,47)=(-(exp(y(47))*params(59)));
g1(52,11)=1;
g1(52,36)=(-(params(60)*exp(y(36))));
g1(53,2)=(-((-((1-params(17)*(1-params(16)))*params(17)*params(20)))/((1-params(17)*(1-params(16)))*y(2)*(1-params(17)*(1-params(16)))*y(2))));
g1(53,22)=1;
g1(54,26)=(-1);
g1(54,54)=1;
g1(55,26)=(-1);
g1(55,55)=1;
g1(56,26)=(-1);
g1(56,56)=1;
g1(57,26)=(-1);
g1(57,57)=1;
g1(58,26)=(-1);
g1(58,58)=1;
g1(59,26)=(-1);
g1(59,59)=1;
g1(60,26)=(-1);
g1(60,60)=1;
g1(61,26)=(-1);
g1(61,61)=1;
g1(62,26)=(-1);
g1(62,62)=1;
g1(63,26)=(-1);
g1(63,63)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end