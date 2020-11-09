function [x1] = GSS(f,xl,xu)
    % Golden Section Search, unconstrained optimization algorithm
    % (Algorithm 11.6.1 of Heer and Maussner, 2009)
    %
    % Input:    f:  pointer to the function f(x)
    %
    %           xl: scalar, lower bound of the interval in which the
    %           maximizer lies
    %
    %           xu: scalar, upper bound of the interval in which the
    %           maximizer lies
    %
    % Output:   x1: the approximate maximizer
    
    % Compute the parameters of the problem
    tol=sqrt(MachEps);
    
    p=(sqrt(5)-1)/2;
    q=1-p;
    
    % Compute the initial interval [a,d] and the points b and c that divide
    % it
    a=xl;
    d=xu;
    b=p*a+q*d;
    c=q*a+p*d;
    
    % Compute the function value at b and c
    fb=f(b);
    fc=f(c);
    
    % Iterate so that the inverval gets smaller
    while abs(d-a)>tol*max([1;(abs(a)+abs(c))])
        if fb<fc % choose [b,d] as next interval
            a=b;
            b=c;
            fb=fc;
            c=p*b+q*d;
            fc=f(c);
        else % choose [a,c] as next interval 
            d=c;
            c=b;
            fc=fb;
            b=p*c + q*a;
            fb=f(b);
        end
    end
 
    if fb>fc
        x1=b;
    else
        x1=c;
    end
    
end
