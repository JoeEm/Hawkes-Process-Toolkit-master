function A = Infectivity_TVHP( T, t, Period, Shift, MaxInfect, Type )

% Generate infectivity matrix A(t) in R^{U*U}
%
% T: the time interval
% t: current time
% Period in R^{U*U}, the predefined period parameter
% Shift in R^{U*U}, the predefined time-shift parameter

switch Type
    case 1
        A = MaxInfect*(1+cos((2*Period./T).*t-Shift));
        if A > ((1/4)*MaxInfect)
            M = 0;
        else
            M = 1;
        end
           A = 0.5*MaxInfect*(1 - (-1).^M);
    case 2
        A = MaxInfect*(1+cos((2*Period./T).*t-Shift));
end

