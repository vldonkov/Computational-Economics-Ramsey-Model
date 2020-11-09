% Definition of utility function

function [c] = rf(z,k1,k2)

    global alpha_coeff delta eta;

    c = z*(k1^alpha_coeff) + (1-delta)*k1 - k2;
    
    if c<0
        c = nan;
    else
        if eta ==1
            c = log(c);
        else
            c=(c^(1-eta))/(1-eta);
        end
    end
        
end