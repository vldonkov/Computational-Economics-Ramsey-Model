function vij = rhs_bellman(z0, x0, x1)
    % Purpose: returns the right-hand side of the Bellman equation
    % rf(z0,x0,x1) + beta_disc*pmat(z0,:)*v1(x1,:)'
    %
    % Input:    z0:  scalar integer, the index of the current shock in the
    %                vector VI_zvec
    %
    %           x0:  scalar, the current value of the endogenous state
    %                variable
    %
    %           x1:  scalar, the next-period value of the endogenous state
    %                variable
    %
    % Output:   vij: the value of the right-hand side of the Bellman
    %                equation at (z0,x0,x1)
    %
    % Remarks: vij is found from linear interpolation (global variable
    %          VI_IP == 1)
    
    global VI_pmat VI_zvec VI_beta_disc VI_xvec VI_ymat VI_IP;
    
    [~,m]=size(VI_pmat);
    vij=rf(VI_zvec(z0),x0,x1);
    if VI_IP==1
        for j=1:m
            vij=vij+VI_beta_disc*VI_pmat(z0,j)*LIP(VI_xvec,VI_ymat(:,j),x1);
        end
    end

end