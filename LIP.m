function [y0] = LIP(xvec,yvec,x0)
    % Purpose: given a function y=f(x) tabulated in xvec and yvec and a
    % point x0, return the function value y0=f(x0) obtained from linear
    % interpolations between x1<x<x2 (as explained on pages 575-576 of Heer
    % and Maussner, 2009)
    %
    % Input:    xvec: n x 1 vector, the tabulated values of the independent
    %                 variable x
    %
    %           yvec: n x 1 vector, the tabulated values of the dependent
    %           variable y
    %
    %           x0:   m x 1 vector, the values of the independent variable
    %           for which y0 is to be computed
    %
    % Output:   y0:   m x 1 vector, see above
    %
    
    n=length(xvec);
    m=length(x0);
    y0=zeros(m,1);
    
    for k=1:m
        if (x0(k)<xvec(1)) || (x0(k)>xvec(n))
            warning("Input out of grid. Procedure will return a missing value. Press any key")
            pause;
            y0 = nan(m,1);
            return;
        end
        if x0(k)==xvec(1)
            y0(k)=yvec(1);
        elseif x0(k)==xvec(n)
            y0(k)=yvec(n);
        else
            j=sum(xvec<=x0(k)); %this determines the lower bracket for x0
            y0(k)=yvec(j)+((yvec(j+1)-yvec(j))/(xvec(j+1)-xvec(j)))*(x0(k)-xvec(j));
        end
    end
    
end