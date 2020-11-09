function [z] = BLIP(xvec,yvec,zmat,x,y)
    % Purpose: Bilinear Interpolation (as explained on pages 576-577 of
    % Heer and Maussner, 2009)
    %
    % Input:    xvec: n x 1 vector, the grid of variable x
    %
    %           yvec: m x 1 vector, the grid of variable y
    %
    %           zmat: n x m matrix, with tabulated function values at
    %           z(i,j) = f(x(i),y(j)), i=1, ..., n, j=1, ..., m
    %
    %           x:    scalar, the x coordinate
    %
    %           y:    scalar, the y coordiante
    %
    % Output:   z:    the interpolated value of f(x,y)
    %
    % Remarks: the elements in xvec and yvec must satisfy xvec(i)<xvec(i+1)
    % for all i and similar for y

    n=length(xvec);
    m=length(yvec);
        
    %first, locate the square that surrounds (x,y)
    if (x<xvec(1)) || (x>xvec(n))
        warning("%f outside of grid! Program stops. Press any key...", x)
        pause;
        z = nan;
        return;
    end
    if (y<yvec(1)) || (y>yvec(m))
        warning("%f outside of grid! Program stops. Press any key...", y)
        pause;
        z = nan;
        return;
    end

    i=sum(xvec<=x);
    j=sum(yvec<=y);
    
    if i==n && j==m
        z=zmat(n,m);
    elseif i==n && j<m
        u=(y-yvec(j))/(yvec(j+1)-yvec(j));
        z=(1.-u)*zmat(n,j)+u*zmat(n,j+1);
    elseif (i<n) && (j==m)
        t=(x-xvec(i))/(xvec(i+1)-xvec(i));
        z=t*zmat(i+1,m)+(1.-t)*zmat(i,m);
    else
        t=(x-xvec(i))/(xvec(i+1)-xvec(i));
        u=(y-yvec(j))/(yvec(j+1)-yvec(j));
        z=(1.-t)*(1.-u)*zmat(i,j)+t*(1.-u)*zmat(i+1,j)+t*u*zmat(i+1,j+1)+(1.-t)*u*zmat(i,j+1);
    end

end