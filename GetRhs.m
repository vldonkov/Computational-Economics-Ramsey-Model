% Right-hand side of Euler equation

function [a] = GetRhs(x)

    global z0 k1 alpha_coeff delta eta;
    
    z1=z0+x; % note that z0=rho*ln(zj)
    z1=exp(z1);
    
    k2=PF(k1,z1);
    c2=z1*(k1^alpha_coeff)+(1-delta)*k1-k2;
    
    a=(c2^(-eta))*(1.-delta+alpha_coeff*z1*(k1^(alpha_coeff-1)));

end