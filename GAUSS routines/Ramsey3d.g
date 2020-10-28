@ ------------------------------------ Ramsey3d.g ----------------------------------

 
  3 February 2008
  9 February 2008 (last changes)
  Alfred Maussner

  Purpose: Solve the stochastic infinite horizon Ramsey model
           of Heer and Maussner, 2nd ed. via iterations over
           the policy function on a discrete grid.

                      1. one-period utility function
                      
                         u(c)=(c^(1-eta)/(1-eta))/(1-eta)
                         
                         or

                         u(C)=ln(c)
                         
                         for eta=1

                      2. technology Y = Z*K^alpha

                      3. transition equation

                         K' = Y + (1-delta)*K - C

                      4. Productivity Shock ln(Z[t]) = rho*ln(Z[t-1])+ eps[t];

  In this version of the program, the previous value function is used
  to initialize the value function for the next, finer grid.

----------------------------------------------------------------------------------- @

new; @ clear memory @

#include SolveDSS.src
#include Tools.src
#include Integration.src

/* Parameters of the model */
alpha=0.27;                   @ elasticity of production with respect to capital    @
delta=0.011;                  @ rate of capital depreciation                        @
beta_disc=0.994;              @ discount factor: 6.5% annual rate of return on capital@
eta=2;                        @ elasticity of marginal utility                    @
rho=0.90;
sigma=0.0072;

/* Parameters of the algorithm */
nz=31;                        @ number of grid points for the productivity shock @
size=5.5;                     @ size of the grid for the productivity shock @
nk=250;                       @ number of grid points for the capital stock @
kmin_g=0.60;                  @ lower bound of the grid for the capital stock @
kmax_g=1.40;                  @ upper bound of the grid for the capital stock @
_VI_IP=0;                     @ =0 without interpolation, =1 linear interpolation, =2 cubic interpolation @
_VI_MPI=1;                    @ =0 without modified policy iteration, =1 with modified policy iteration @
_VI_MPI_K=30;                 @ number of iterations over the policy function @
_VI_nc=50;                    @ number of iterations with unchanged policy function before iterations are terminated @
_cumulative=1;                @ in successive computations, compute total run time, =0 do not @

/* Parameters for the computation of Euler equation residuals */
kmin_e=0.8;               @ kmin_e*kstar is the lower bound @
kmax_e=1.2;               @ kmax_e*kstar is the upper bound @
nobs_e=200;               @ the number of residuals to be computed @


/* Open file and write initial information */
output file=Ramsey3d_2lam1.txt reset;
MyDate;
?"Parameters of the model:";
?"alpha = " ftos(alpha,"*.*lf",8,4);
?"beta_disc  = " ftos(beta_disc, "*.*lf",8,4);
?"delta = " ftos(delta,"*.*lf",8,4);
?"eta   = " ftos(eta,  "*.*lf",8,4);
?"rho   = " ftos(rho,  "*.*lf",8,4);
?"sigma = " ftos(sigma,"*.*lf",8,4);
?"";
?"Parameters of the algorithm:";
?"kmin_g= " ftos(kmin_g, "*.*lf",8,3);
?"kmax_g= " ftos(kmax_g, "*.*lf",8,3);
?"kmin_e= " ftos(kmin_e, "*.*lf",8,3);
?"kmax_e= " ftos(kmax_e, "*.*lf",8,3);
?"nobs_e= " ftos(nobs_e,"*.*lf",8,0);
?"nz    = " ftos(nz    ,"*.*lf",8,0);
?"";
?"eps   = " ftos(_VI_eps ,"*.*lf",7,4);
?"stop  = " ftos(_VI_nc  ,"*.*lf",7,0);
?"IP    = " ftos(_VI_IP ,"*.*lf",7,0);
?"MPI   = " ftos(_VI_MPI,"*.*lf",7,0);
if _VI_MPI;
?"K     = " ftos(_VI_MPI_K,"*.*lf",7,0);
endif;

output off;

/* Compute Markov chain approximation */
{zgrid,pmat}=MarkovAR(size,nz,rho,sigma);

zgrid=exp(zgrid);

zmin=zgrid[1];
zmax=zgrid[nz];

zmin_i=rho*ln(zmin)+sqrt(2)*sigma*(-1.650680123);
zmax_i=rho*ln(zmax)+sqrt(2)*sigma*(1.650680123);


kmin=((1-beta_disc*(1-delta))/(alpha*beta_disc*zmin))^(1/(alpha-1));
kmax=((1-beta_disc*(1-delta))/(alpha*beta_disc*zmax))^(1/(alpha-1));

/* Compute stationary solution of deterministic model and intialize the value function */
kstar=((1- beta_disc*(1-delta))/alpha*beta_disc)^(1/(alpha-1));  @ stationary capital stock                          @
cstar=kstar^alpha - delta*kstar;                       @ stationary level of consumption                   @


kmin_g=kmin_g*kmin;
kmax_g=kmax_g*kmax;

kmin_e=kmin_e*kstar;
kmax_e=kmax_e*kstar;

/* Vector with different values of nk */
nvec={250, 1000, 10000};

/* Iterations of nvec start here */
lmax=rows(nvec);

policy=arrayinit(lmax|nobs_e|nobs_e,0);
  emat=arrayinit(lmax|nobs_e|nobs_e,0);

/* Different ways to initialize v0, uncomment to try the others */
nk=nvec[1];
v0=rf(1,kstar,kstar)/(1-beta_disc);  @ stationary solution @
v0=ones(nk,nz).*v0;

@     v0=zeros(nk,nz); @        @ the zeros function  @
@    for i (1,nk,1); for j (1,nz,1); knext=kgrid[i]; v0[i,j]=rf(zgrid[j],kgrid[i],knext); endfor; endfor;@ @ maintain the given capital stock @

/* Iterations over different nk start here */
s2=0; @ stores time needed to obtain initial v from previous v @

tottime=zeros(lmax,1);

for l (1,lmax,1);
    nk=nvec[l];
    kgrid=seqa(kmin_g,(kmax_g-kmin_g)/(nk-1),nk);
    
    @ outcommend if necessary @
    @ if l==1; _VI_IP=0; else; _VI_IP=2; endif; @

    /* Solve for the policy function */
    s1=hsec;
    {v1,hmati}=SolveVIS(beta_disc,kgrid,zgrid,pmat,v0);
    s1=hsec-s1;
    
    if _VI_IP/=0;
     hmat=hmati;
    else;
        hmat=zeros(nk,nz);
        for i (1,nk,1);
            for j (1,nz,1);
                hmat[i,j]=kgrid[hmati[i,j]];
            endfor;
        endfor;
    endif;

    /* Computation of Euler equation residuals */
    kvec=seqa(kmin_e,(kmax_e-kmin_e)/(nobs_e-1),nobs_e);
    zvec=seqa(0.95,0.1/(nobs_e-1),nobs_e);
    clear z0, k1;
    eer=Euler(kvec,zvec);
    emat[l,1:nobs_e,1:nobs_e]=eer;
    emax=maxc(maxc(abs(eer)));
    tottime[l]=s1+s2;
    test="nk= " $+ ftos(nk,"*.*lf",10,0) $+  " Run time= " $+ etstr(s1+s2);
    if _cumulative; test=test $+ "Cumulative run time= " $+etstr(sumc(tottime[1:l])); endif;
    test=test $+ " EER= " $+ ftos(emax,"*.*lf",10,8);
    output on;
    @?"_VI_IP= " ftos(_VI_IP,"*.*lf",2,0);@
    ?test;
    output off;
    /* computation of policy function */
    for i (1,nobs_e,1);
        for j (1,nobs_e,1);
            policy[l,i,j]=BLIP(kgrid,zgrid,hmat,kvec[i],zvec[j]);
        endfor;
    endfor;
    /* New initial v0 */    
    if l<lmax;
        if nk==nvec[l+1];
            v0=v1;
        else;
            locate 15,5; ?"Compute new initial v0";
            s2=hsec;
            nk1=nvec[l+1];
            kgnew=seqa(kmin_g,(kmax_g-kmin_g)/(nk1-1),nk1);   
            kgnew[nk1]=kgrid[nk];kgnew[1]=kgrid[1]; 
            v0=zeros(nk1,nz);    
            for j (1,nz,1);
                locate 15,30;?ftos(j,"*.*lf",9,0);
                v0[.,j]=LIP(kgrid,v1[.,j],kgnew);
            endfor;
            s2=hsec-s2;
        endif;
    endif;
endfor;
save emat, policy, kvec, zvec;
ende:
end;
   

/* ----------------------------------- Subroutines ----------------------- */
   
proc(1)=rf(z,k1,k2);  @ defines the utility function @

  local c;
  c=z*(k1^alpha) + (1-delta)*k1 - k2;

  if c<0.0;
     c=miss(1,1);
  else;
     if eta==1;
       c=ln(c);
     else;
       c=(c^(1-eta))/(1-eta);
     endif;
  endif;

retp(c);

endp;

/* Policy function via linear interpolation: kgrid, zgrind, and hmat must be given as a global matrices */

proc(1)=PF(k,z);

    local knext;

    knext=BLIP(kgrid,zgrid,hmat,k,z);

    retp(knext);

endp;


/* Euler equation residuals */
proc(1)=Euler(kvec,zvec);

    local n, m, eer, i, j, c0, c1, rhs;
    external proc GH_INT4, GetRhs;
    external matrix z0, k1;

     n=rows(kvec);
     m=rows(zvec);
   eer=zeros(n,m);


    for i (1,n,1);
        for j (1,m,1);
			z0=rho*ln(zvec[j]);
            k1=PF(kvec[i],zvec[j]);
            c0=zvec[j]*(kvec[i]^alpha)+(1-delta)*kvec[i]-k1;            
            rhs=GH_INT4(&GetRhs,sigma);
			rhs=rhs*beta_disc;
            c1=(rhs^(-1/eta));
			eer[i,j]=(c1/c0)-1;
        endfor;
    endfor;

    retp(eer);

	endp;

/* Rhs of Eulerequation */
proc(1)=GetRhs(x);

 local z1, c2, k2;
 external matrix z0, k1;
	
	z1=z0+x;  @ note that z0=rho*ln(zj) @    
    z1=exp(z1);

    k2=PF(k1,z1);
    c2=z1*(k1^alpha)+(1-delta)*k1-k2;
    
    retp((c2^(-eta))*(1.-delta+alpha*z1*(k1^(alpha-1))));

endp;

