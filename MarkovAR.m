function [zt,p] = MarkovAR(size_z,m,rho,sigma)
    % Purpose: Approximate AR(1)-Process by Markov chain (Algorithm 12.2.1
    % of Heer and Maussner, 2009)
    %
    % Input:    size_z:  scalar, the mulitple of the unconditional standard
    %                  deviation of the AR(1) process used to define the
    %                  grid size
    %
    %           m:     integer scalar, the number of grid points
    %
    %           rho:   scalar, the autoregressive parameter of the process
    %
    %           sigma: scalar, the standard deviation of the innovations
    %
    %
    % Output:   zt:     m x 1 vector, the grid approximating the process
    %
    %           p:      m x m matrix of transition probabilities

    sigmaz=sqrt(sigma^2/(1-rho^2));
    zbar=size_z*sigmaz;
    N=m;
    start_val=-zbar;
    inc = 2*zbar/(m-1);
    stop_val = (N-1)*inc + start_val;
    zt=start_val:inc:stop_val;
    
    p=zeros(m);
    i=1;
    while i<=m
        p(i,1)=normcdf((zt(1)-rho*zt(i)+(zt(2)-zt(1))/2)/sigma,0,1);
        
        j=2;
        while j<=(m-1)
            p(i,j)=normcdf((zt(j)-rho*zt(i)+(zt(j)-zt(j-1))/2)/sigma,0,1)-...
                normcdf((zt(j)-rho*zt(i)-(zt(j)-zt(j-1))/2)/sigma,0,1);
            j=j+1;
        end
        p(i,m)=1-sum(p(i,1:(m-1)));
        i=i+1;
    end
    
end