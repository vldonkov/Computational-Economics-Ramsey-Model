% Clear workspace, command window and close all figures 
clear all;
close all;
clc;

                % Stochastic Infinite-Horizon Ramsey Model

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
%            Y = Z*K^alpha_coeff
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
global alpha_coeff delta eta beta_disc rho sigma;
alpha_coeff=0.27; % elasticity of production with respect to capital
delta=0.011; % rate of capital depreciation
beta_disc=0.994; % discount factor: 6.5% annual rate of return on capital
eta=2; % elasticity of marginal utility
rho=0.90; % parameter of continious-valued AR(1)-process
sigma=0.0072; % parameter of continious-valued AR(1)-process

% Parameters of the algorithm
nz=9; % number of grid points for the productivity shock
size_z_gr=5.5; % size of the grid for the productivity shock
nvec=[5,25,250]; % Vector with different values of the number of
                       % grid points for the capital stock
kmin_g=0.60; % lower bound of the grid for the capital stock
kmax_g=1.40; % upper bound of the grid for the capital stock
global VI_IP VI_nc;
VI_IP=1; % =0 without interpolation, =1 linear interpolation
VI_nc=50; % stop if number of consecutive iterations with unchanged 
          % indices of policy function exceeds this number
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

% Declaration of some global variables used throughout subroutines
global z0 k1 kgrid zgrid hmat;


% Open file and write initial information
file = fopen("RamseyModel.txt","w");
fprintf(file,"%s\n",datetime('now'));
fprintf(file,"\n");
fprintf(file,"Parameters of the model:\n");
fprintf(file,"alpha_coeff = %.4f\n",alpha_coeff);
fprintf(file,"beta_disc   = %.4f\n",beta_disc);
fprintf(file,"delta       = %.4f\n",delta);
fprintf(file,"eta         = %.4f\n",eta);
fprintf(file,"rho         = %.4f\n",rho);
fprintf(file,"sigma       = %.4f\n",sigma);
fprintf(file,"\n");
fprintf(file,"Parameters of the algorithm:\n");
fprintf(file,"nvec      = [%d, %d, %d]\n",nvec);
fprintf(file,"kmin_g    = %.3f\n",kmin_g);
fprintf(file,"kmax_g    = %.3f\n",kmax_g);
fprintf(file,"kmin_e    = %.3f\n",kmin_e);
fprintf(file,"kmax_e    = %.3f\n",kmax_e);
fprintf(file,"nobs_e    = %d\n",nobs_e);
fprintf(file,"nz        = %d\n",nz);
fprintf(file,"size_z_gr = %.3f\n",size_z_gr);
fprintf(file,"\n");
fprintf(file,"Important parameters of function SolveVIS:\n");
fprintf(file,"Max  = %d\n",VI_Max);
fprintf(file,"eps  = %.4f\n",VI_eps);
fprintf(file,"stop = %d\n",VI_nc);
fprintf(file,"IP   = %d\n",VI_IP);
fprintf(file,"\n");
fclose(file);


% Here beings the actual algorithm of the stochastic Ramsey model

% Compute Markov chain approximation
[zgrid,pmat]=MarkovAR(size_z_gr,nz,rho,sigma);

zgrid=exp(zgrid);

zmin=zgrid(1);
zmax=zgrid(nz);

zmin_i=rho*log(zmin)+sqrt(2)*sigma*(-1.650680123);
zmax_i=rho*log(zmax)+sqrt(2)*sigma*(1.650680123);


kmin=((1-beta_disc*(1-delta))/(alpha_coeff*beta_disc*zmin))^(1/(alpha_coeff-1));
kmax=((1-beta_disc*(1-delta))/(alpha_coeff*beta_disc*zmax))^(1/(alpha_coeff-1));

% Compute stationary solution of deterministic model and intialize the
% value function
kstar=((1- beta_disc*(1-delta))/alpha_coeff*beta_disc)^(1/(alpha_coeff-1)); % Stationary capital stock
cstar=kstar^alpha_coeff - delta*kstar; % Stationary level of consumption


kmin_g=kmin_g*kmin;
kmax_g=kmax_g*kmax;

kmin_e=kmin_e*kstar;
kmax_e=kmax_e*kstar;

% Different ways to initialize v0 (uncomment the one you want to try)

% Option 1: zero matrix
% nk=nvec(1);
% v0=zeros(nk,nz);

% Option 2: maintain the given capital stock
% nk=nvec(1);
% N=nk;
% start_val=kmin_g;
% inc = (kmax_g-kmin_g)/(nk-1);
% stop_val = (N-1)*inc + start_val;
% kgrid=start_val:inc:stop_val;
% v0=zeros(nk,nz);
% for i=1:nk
%     for j=1:nz
%         knext=kgrid(i);
%         v0(i,j)=rf(zgrid(j),kgrid(i),knext);
%     end
% end

% Option 3: stationary solution of the deterministic growth model
nk=nvec(1);
v0=rf(1,kstar,kstar)/(1-beta_disc);
v0=ones(nk,nz).*v0;


% Iterations of nvec start here
lmax=length(nvec);

policy = zeros(nobs_e,nobs_e,lmax);
emat = zeros(nobs_e,nobs_e,lmax);

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
    
    % Computation of Euler equation residuals
    N=nobs_e;
    start_val=kmin_e;
    inc=(kmax_e-kmin_e)/(nobs_e-1);
    stop_val=(N-1)*inc + start_val;
    kvec=start_val:inc:stop_val; 
    start_val=0.95;
    inc=0.1/(nobs_e-1);
    stop_val=(N-1)*inc + start_val;
    zvec=start_val:inc:stop_val;
    z0=0;
    k1=0;
    eer=Euler(kvec,zvec);
    emat(1:nobs_e,1:nobs_e,l)=eer;
    emax=max(max(abs(eer)));
    tottime(l)=s1+s2;
    file = fopen("RamseyModel.txt","a+");
    fprintf(file,"nk = %d\n",nk);
    fprintf(file,"Run time = %d minutes and %.2f seconds\n",floor((s1+s2)/60),rem((s1+s2),60));
    if cumulative
        fprintf(file,"Cumulative run time = %d minutes and %.2f seconds\n",floor(sum(tottime(1:l))/60),rem(sum(tottime(1:l)),60));
    end
    fprintf(file,"EER = %e\n",emax);
    fprintf(file,"\n");
    fclose(file);
    % computation of policy function
    for i=1:nobs_e
        for j=1:nobs_e
            policy(i,j,l)=BLIP(kgrid,zgrid,hmat,kvec(i),zvec(j));
        end
    end
    % New initial v0
    if l<lmax
        if nk==nvec(l+1)
            v0=v1;
        else
            fprintf("\n")
            tic;
            nk1=nvec(l+1);
            N=nk1;
            start_val=kmin_g;
            inc=(kmax_g-kmin_g)/(nk1-1);
            stop_val=(N-1)*inc + start_val;
            kgnew=start_val:inc:stop_val; 
            kgnew(nk1)=kgrid(nk);
            kgnew(1)=kgrid(1);
            v0=zeros(nk1,nz);
            for j=1:nz
                fprintf("Compute new initial v0        %d\n",j)
                v0(:,j)=LIP(kgrid,v1(:,j),kgnew);
            end
            s2=toc;
        end
    end 
end
save('emat.mat','emat');
save('policy.mat','policy');
save('kvec.mat','kvec');
save('zvec.mat','zvec');