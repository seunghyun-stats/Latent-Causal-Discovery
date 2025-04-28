function st = ST(theta, gamma)
% this function is used to calculate ST function in updating d
%
% @param theta: item parameter matrix
% @param gamma: tuning parameter for TLP
%
% @return st: value of ST function

    if abs(theta) > gamma
        st = (abs(theta) - gamma) * sign(theta);
    else
        st = 0;
    end

end