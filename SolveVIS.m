function [v1,xz]=SolveVIS(beta_disc,xvec,zvec,pmat,v0)
    % Purpose:   Computes the policy function for a stochastic Dynamic
    %            General Equilibrium model with one endogenous state
    %            variable x and one exogenous shock z.
    %
    % Algorithm: The policy function is computed on a matrirx of of n by m
    %            points, where n is the number of points in the grid for x
    %            and m the number of points in the grid for z.
    %
    %            The method employed is iteration over the value function.
    %            The user can choose whether linear interpolation between
    %            the grid points of x should be used. For linear
    %            interpolation the global variable VI_IP must be set to 1.
    %            The default is 0, i.e., no interpolation is used. In this
    %            case, the solution matrix xz is a matrix of intergers
    %            (i,j), where the (i,j) element of this matrix points to
    %            the index of the x that gives the optimal next-period
    %            value of the endogenous state variable.
    %
    %            The  algorithm makes use of the concavity of the value
    %            function and the montonicity of the policy function.
    %            Specifically, a binary search algorithm is being used in
    %            order to locate the maximum on the right-hand side of the
    %            Bellman equation.
    %
    %
    % Input:     beta_disc: the discount factor from the agent's life-time
    %                       utility function
    %
    %            xvec:      n x 1 column vector with elements in ascending
    %                       order that represent valid elements of the
    %                       variable x
    %
    %            zvec:      m x 1 column vector that represents the
    %                       elements of the Markov chain, which represents
    %                       the stochastic variable z
    %
    %            pmat:      m x m matrix, the transition matrix of the
    %                       Markov chain
    %
    %            v0:        n x m matrix with the elements of the initial
    %                       value function
    %
    %
    % Output:    v1:        n x m matrix, the value function at the
    %                       solution
    %
    %            xz:        n x m matrix, if i=1, ..., n is the index of a
    %                       point from xvec and j=1, ..., m is the index of
    %                       a point from zvec, xz(i,j)=jstar is the index
    %                       of xvec that gives the optimal next-period
    %                       value of the endogeneous state variable x.
    %
    % Remarks:   The one-period utility function rf(z,x,x') must be
    %            programmed in a procedure lhs=rf(z,x,x') with name rf and
    %            the specified ordering of the variables. lhs returns the
    %            one-period utility, if the current state of the system is
    %            given by (x,z) and if a the next-period state x' is
    %            chosen.
    
    global VI_IP VI_nc;
    global VI_xvec VI_ymat VI_zvec VI_pmat VI_beta_disc VI_xex VI_zex;
    global VI_eps VI_Max;
    
    % Step 1: Initialize
    eps1=VI_eps*(1-beta_disc); % convergence criteria
    
    VI_beta_disc=beta_disc;

    nx=length(xvec); % number of grid points in xvec
    nz=length(zvec); % number of grid points in zvec
    
    if VI_IP==0
        h1=ones(nx,nz); % intial policy function
    end

    h2=zeros(nx,nz); % new policy function
    w =zeros(3,1);
    
    v1=v0; % old policy function
    v2=zeros(nx,nz); % new policy function
    dv=1;
    nc=0;
    
    if VI_IP==1
        VI_xvec=xvec;
        VI_zvec=zvec;
        VI_ymat=v0;
        VI_pmat=pmat;
        VI_beta_disc=beta_disc;
    end
    
    % Step 2: Iterate over the value function
    t=1;
    while (t<=VI_Max) && (dv>=eps1) && (nc<=VI_nc) % begin loop over value function
        
        for j=1:nz % begin loop over zvec
            if VI_IP==1
                VI_zex=j;
            end
            js=1;
            for i=1:nx % begin loop over xvec
                if VI_IP==1
                    VI_xex=xvec(i);
                end 
                jmin=js;
                jmax=nx;
                while (jmax-jmin)>2 % the next lines implement a binary
                                    % search algorithm
                    jl=floor((jmin+jmax)/2);
                    ju=jl+1;
                    w(1)=rf(zvec(j),xvec(i),xvec(jl))+beta_disc*(pmat(j,:)*(v1(jl,:)'));
                    w(2)=rf(zvec(j),xvec(i),xvec(ju))+beta_disc*(pmat(j,:)*(v1(ju,:)'));
                    if w(2)>w(1)
                        jmin=jl;
                    else
                        jmax=ju;
                    end
                end
                w(1)=rf(zvec(j),xvec(i),xvec(jmin))+beta_disc*(pmat(j,:)*(v1(jmin,:)'));
                if jmax>jmin
                    w(2)=rf(zvec(j),xvec(i),xvec(jmin+1))+beta_disc*(pmat(j,:)*(v1(jmin+1,:)'));
                else
                    w(2)=w(1);
                end
                w(3)=rf(zvec(j),xvec(i),xvec(jmax))+beta_disc*(pmat(j,:)*(v1(jmax,:)'));
                [~,js]=max(w);
                if VI_IP==0
                    v2(i,j)=w(js);
                end
                js=jmin+js-1;
                
                % The next lines implement linear interpolation between grid points
                if VI_IP==1
                    if js==1 % boundary optimum, ax=bx=a(1)
                        ax=xvec(1);
                        bx=ax+eps1*(xvec(2)-xvec(1));
                        if rhs_bellman(j,xvec(i),bx)<rhs_bellman(j,xvec(i),ax)
                            h2(i,j)=xvec(1);
                        else
                            h2(i,j)=GSS(@VI_valuefunction,xvec(1),xvec(2));
                        end
                    elseif js==nx % boundary optimum, bx=cx=a(n)
                        cx=xvec(nx);
                        bx=cx-eps1*(xvec(nx)-xvec(nx-1));
                        if rhs_bellman(j,xvec(i),bx)<rhs_bellman(j,xvec(i),cx)
                            h2(i,j)=xvec(nx);
                        else
                            h2(i,j)=GSS(@VI_valuefunction,xvec(nx-1),xvec(nx));
                        end
                    else
                        h2(i,j)=GSS(@VI_valuefunction,xvec(js-1),xvec(js+1));
                    end
                    
                    v2(i,j)=rhs_bellman(j,xvec(i),h2(i,j));
                else
                    h2(i,j)=js;
                end
                
            end % end loop over xvec    
        end % end loop over zvec
        
        if VI_IP==0
            % compute stopping criterium 2
            di=sum(sum(h2 ~= h1));
            if di>=1
                nc=0;
            else
                nc=nc+1;
            end
            h1=h2;
        end
        dv=max(max(abs(v2-v1)));
        
        clc;
        fprintf("Iteration #= %d\n", t)
        fprintf("Largest element in v1-v0= %f\n", dv)
        if VI_IP==0
            fprintf("# of indices that have changed: %d\n", di)
            fprintf("# of consecutive iterations with constant policy function= %d\n", nc)
        end
        v1=v2;
        if VI_IP==1
            VI_ymat=v1;
        end
        t=t+1;    
    end
    fprintf("\n")
    if t>VI_Max
        warning("Maximum number of iterations exceeded. Change VI_Max!")
        warning("The computed solution may be inaccurate.Press any key...")
        pause;
    end
    if VI_IP==0
        if min(min(h1))==1
            warning("Policy function hits lower bound of grid")
        end
        if max(max(h1))==nx
            warning("Policy function hits upper bound of grid")
        end
    else
        if min(min(h2))==xvec(1)
            warning("Policy function hits lower bound of grid")
        end
        if max(max(h2))==xvec(nx)
            warning("Policy function hits upper bound of grid")
        end
    end
    
    xz=h2;
    
end