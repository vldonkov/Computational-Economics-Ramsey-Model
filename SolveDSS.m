
%% Clear workspace, command window and close all figures 
clear all;
close all;
clc;


                %% Stochastic Infinite-Horizon Ramsey Model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Purpose: Solve a stochastic infinite-horizon Ramsey model via value
%            function iterations with/without linear interpolation on a
%            discrete grid.
%
%            1. one-period utility function
%                   
%            u(c)=(c^(1-eta)-1)/(1-eta)
%
%            or
%
%            u(c)=ln(c)
%                     
%            for eta=1
%
%            2. technology
%
%            Y = Z*K^alpha
%
%            3. transition equation
%
%            K' = Y + (1-delta)*K - C
%
%            4. productivity shock
%
%            ln(Z(t)) = rho*ln(Z(t-1))+ sigma*eps(t)
%
%   Remarks: In this version of the program, the previous value function is
%            used to initialize the value function for the next, finer
%            grid.


% Parameters of the model
alpha=0.27; % elasticity of production with respect to capital
delta=0.011; % rate of capital depreciation
beta_disc=0.994; % discount factor: 6.5% annual rate of return on capital
eta=2; % elasticity of marginal utility
rho=0.90; % parameter of continious-valued AR(1)-process
sigma=0.0072; % parameter of continious-valued AR(1)-process

% Parameters of the algorithm
nz=31; % number of grid points for the productivity shock
size=5.5; % size of the grid for the productivity shock
nvec=[250 1000 10000]; % Vector with different values of the number of
                       % grid points for the capital stock
kmin_g=0.60; % lower bound of the grid for the capital stock
kmax_g=1.40; % upper bound of the grid for the capital stock
global VI_IP VI_nc;
VI_IP=0; % =0 without interpolation, =1 linear interpolation
VI_nc=50; % number of iterations with unchanged policy function before
          % iterations are terminated
cumulative=1; % in successive computations, compute total run time,
              % =0 do not

% Parameters for the computation of Euler equation residuals
kmin_e=0.8; % kmin_e*kstar is the lower bound
kmax_e=1.2; % kmax_e*kstar is the upper bound
nobs_e=200; % the number of residuals to be computed

% The next lines set default values for the global variables used by the
% function SolveVIS
global VI_xvec VI_ymat VI_zvec VI_pmat VI_beta_disc VI_xex VI_zex;
global VI_eps VI_Max;
VI_xvec=0; % stores the x values used in interpolation
VI_ymat=0; % stores the y values uses in interpolation
VI_zvec=0; % stores the z values used in interpolation
VI_pmat=0; % stores the transition matrix used in interpolation
VI_beta_disc=0; % stores information used for linear interpolation
VI_xex=0; % stores information to write the value function as a function of
          % x(i) alone
VI_zex=0; % stores information to write the value function as a function of
          % x(i) alone
VI_eps=0.01; % stopping criterium
VI_Max=1000; % maximal number of iterations


% Open file and write initial information
file = fopen("RamseyModel.txt","w");
fprintf(file,"%s\n",datetime('now'));
fprintf(file,"\n");
fprintf(file,"Parameters of the model:\n");
fprintf(file,"alpha     = %.4f\n",alpha);
fprintf(file,"beta_disc = %.4f\n",beta_disc);
fprintf(file,"delta     = %.4f\n",delta);
fprintf(file,"eta       = %.4f\n",eta);
fprintf(file,"rho       = %.4f\n",rho);
fprintf(file,"sigma     = %.4f\n",sigma);
fprintf(file,"\n");
fprintf(file,"Parameters of the algorithm:\n");
fprintf(file,"nvec   = [%d, %d, %d]\n",nvec);
fprintf(file,"kmin_g = %.3f\n",kmin_g);
fprintf(file,"kmax_g = %.3f\n",kmax_g);
fprintf(file,"kmin_e = %.3f\n",kmin_e);
fprintf(file,"kmax_e = %.3f\n",kmax_e);
fprintf(file,"nobs_e = %d\n",nobs_e);
fprintf(file,"nz     = %d\n",nz);
fprintf(file,"size   = %.3f\n",size);
fprintf(file,"\n");
fprintf(file,"Important parameters of function SolveVIS:\n");
fprintf(file,"Max  = %d\n",VI_Max);
fprintf(file,"eps  = %.4f\n",VI_eps);
fprintf(file,"stop = %d\n",VI_nc);
fprintf(file,"IP   = %d\n",VI_IP);
fclose(file);


% Here beings the actual algorithm of the stochastic Ramsey model

% Compute Markov chain approximation
[zgrid,pmat]=MarkovAR(size,nz,rho,sigma);

zgrid=exp(zgrid);

zmin=zgrid(1);
zmax=zgrid(nz);

zmin_i=rho*log(zmin)+sqrt(2)*sigma*(-1.650680123);
zmax_i=rho*log(zmax)+sqrt(2)*sigma*(1.650680123);


kmin=((1-beta_disc*(1-delta))/(alpha*beta_disc*zmin))^(1/(alpha-1));
kmax=((1-beta_disc*(1-delta))/(alpha*beta_disc*zmax))^(1/(alpha-1));

% Compute stationary solution of deterministic model and intialize the
% value function
kstar=((1- beta_disc*(1-delta))/alpha*beta_disc)^(1/(alpha-1)); % Stationary capital stock
cstar=kstar^alpha - delta*kstar; % Stationary level of consumption


kmin_g=kmin_g*kmin;
kmax_g=kmax_g*kmax;

kmin_e=kmin_e*kstar;
kmax_e=kmax_e*kstar;

% Iterations of nvec start here
lmax=length(nvec);

policy = zeros(nobs_e,nobs_e,lmax);
emat = zeros(nobs_e,nobs_e,lmax);

% Different ways to initialize v0, uncomment to try the others
nk=nvec(1);
v0=rf(1,kstar,kstar)/(1-beta_disc); % stationary solution
v0=ones(nk,nz).*v0;

%v0=zeros(nk,nz); % the zeros function
%for i=1:nk % maintain the given capital stock
%    for j=1:nz
%        knext=kgrid(i);
%        v0(i,j)=rf(zgrid(j),kgrid(i),knext);
%    end
%end

% Iterations over different nk start here
s2=0; % stores time needed to obtain initial v from previous v

tottime=zeros(lmax,1);

for l=1:lmax
    nk=nvec(l);
    N=nk;
    start_val=kmin_g;
    inc = (kmax_g-kmin_g)/(nk-1);
    stop_val = (N-1)*inc + start_val;
    kgrid=start_val:inc:stop_val;
    
    % Solve for the policy function
    tic;
    [v1,hmati]=SolveVIS(beta_disc,kgrid,zgrid,pmat,v0);
    s1=toc;
    
    if VI_IP==1
        hmat=hmati;
    else
        hmat=zeros(nk,nz);
        for i=1:nk
            for j=1:nz
                hmat(i,j)=kgrid(hmati(i,j));
            end
        end
    end
    
    %Continue here
    
    
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                                %% SolveVIS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    %                       a point from zvec, xz[i,j]=jstar is the index
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
                        cx=xvec(2);
                        if rhs_bellman(j,xvec(i),bx)<rhs_bellman(j,xvec(i),ax)
                            h2(i,j)=xvec(1);
                        else
                            h2(i,j)=GSS(@VI_valuefunction,xvec(1),xvec(2));
                        end
                    elseif js==nx % boundary optimum, bx=cx=a(n)
                        ax=xvec(nx-1);
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
            printf("# of indices that have changed: %d\n", di)
            printf("# of consecutive iterations with constant policy function= %d\n", nc)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                        %Define the utility function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c] = rf(z,k1,k2)

    c = z*(k1^alpha) + (1-delta)*k1 - k2;
    
    if c<0.0
        c = missing;
    else
        if eta ==1
            c = log(c);
        else
            c=(c^(1-eta))/(1-eta);
        end
    end
        
end





    