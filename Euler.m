% Euler equation residuals

function [eer] = Euler(kvec,zvec)

    global z0 k1 alpha_coeff delta beta_disc eta rho sigma;
    
    n=length(kvec);
    m=length(zvec);
    eer=zeros(n,m);
    
    for i=1:n
        for j=1:m
            z0=rho*log(zvec(j));
            k1=PF(kvec(i),zvec(j));
            c0=zvec(j)*(kvec(i)^alpha_coeff)+(1-delta)*kvec(i)-k1;
            rhs=GH_INT4(@GetRhs,sigma);
            rhs=rhs*beta_disc;
            c1=(rhs^(-1/eta));
            eer(i,j)=(c1/c0)-1;
        end
    end

end