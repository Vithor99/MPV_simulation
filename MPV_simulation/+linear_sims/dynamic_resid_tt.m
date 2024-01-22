function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 15);

T(1) = (params(44)/params(45))^params(7);
T(2) = params(41)^(-params(6));
T(3) = params(40)*y(42)+T(2)*(y(43)*params(6)-y(56));
T(4) = params(25)*(params(4)*params(45)*(1-params(16))*params(1)-params(18)*params(44)*(params(48)+params(16)*(1-params(48))*params(22))/params(49));
T(5) = params(44)*params(25)*(1-params(16)*params(22)+params(16)*params(22)/params(48));
T(6) = params(11)/(1-params(11))*(1-params(17)*(1-params(11)))/(1+params(21)/params(5)-params(21));
T(7) = (params(44)/params(45))^(1+params(7));
T(8) = (params(48)+params(16)*(1-params(48))*params(22))*T(7);
T(9) = params(48)*((1-params(16)*params(22)-(1-params(16))*params(1))/(params(48)+(1-params(48))*(params(16)*params(22)+(1-params(16))*params(1)))-((1-params(16)*params(22))*T(7)-(1-params(16))*params(1))/((1-params(16))*(1-params(48))*params(1)+T(8)));
T(10) = (params(41)*params(43)/params(8))^params(10);
T(11) = params(2)^(1+params(10));
T(12) = params(19)^(1+params(10));
T(13) = params(40)^(1/params(5));
T(14) = params(44)*params(16)*params(22)*params(25)/params(51);
T(15) = params(48)*(1-(1-params(16))*params(1)-params(16)*params(22))/(params(48)+(1-params(48))*(params(16)*params(22)+(1-params(16))*params(1)));

end
