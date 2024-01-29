function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 15);

T(1) = params(11)/(1-params(11))*(1-(1-params(11))*params(17))/(1+params(20)/params(5)-params(20));
T(2) = (params(42)/params(43))^params(7);
T(3) = params(39)^(-params(6));
T(4) = params(38)*y(22)+T(3)*(y(23)*params(6)-y(38));
T(5) = params(24)*(params(4)*params(43)*(1-params(16))*params(1)-params(18)*params(42)*(params(46)+params(16)*(1-params(46))*params(21))/params(47));
T(6) = params(42)*params(24)*(1-params(16)*params(21)+params(16)*params(21)/params(46));
T(7) = (params(42)/params(43))^(1+params(7));
T(8) = (params(46)+params(16)*(1-params(46))*params(21))*T(7);
T(9) = params(46)*((1-params(16)*params(21)-(1-params(16))*params(1))/(params(46)+(1-params(46))*(params(16)*params(21)+(1-params(16))*params(1)))-((1-params(16)*params(21))*T(7)-(1-params(16))*params(1))/((1-params(16))*(1-params(46))*params(1)+T(8)));
T(10) = (params(39)*params(41)/params(8))^params(10);
T(11) = params(2)^(1+params(10));
T(12) = params(3)^(1+params(10));
T(13) = params(38)^(1/params(5));
T(14) = params(42)*params(16)*params(21)*params(24)/params(49);
T(15) = params(46)*(1-(1-params(16))*params(1)-params(16)*params(21))/(params(46)+(1-params(46))*(params(16)*params(21)+(1-params(16))*params(1)));

end
