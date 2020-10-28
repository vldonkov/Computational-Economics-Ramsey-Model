@ ------------------------------------------- Ramsey2d.g ------------------------------------

   Alfred Maußner

   16 January 2008
   25 January 2008 (cubic spline interpolation added)
   26 January 2008 (modified policy iterations with sparse matrix routines)
   
   Test various ways to obtain the solution of the infinite horizon Ramsey model
   of Heer and Maussner, 2nd, ed. via value function iteration.

   

----------------------------------------------------------------------------------------------------- @

new; 

/* Structure to store information and parameter settings for the algorithms */
struct DummyS
{string  MethodName; 
 scalar  VI_PI;     @ =1 full policy function iteration, =2 partial policy function iteration @
 scalar  VI_PI_K;   @ number of iterations in partial policy function iteration @
 scalar  IP;        @ =0 without interpolation, =1: linear interpolation, =2 cubic spline @ 
 scalar  BS;        @ =1: binary search, =0: sequential search @
 scalar  eps;       @ convergence between successive value functions @
 scalar  stop;      @ stop after stop iterations with unchanged indices of policy function @
 scalar  maxit;     @ maximum number of iterations @
};

struct DummyS Method;
Method=reshape(Method,6,1);
Method[1].MethodName="Simple value function iteration";
Method[1].VI_PI=0;
Method[1].IP=0;
Method[1].BS=1;
Method[1].eps=0.01;
Method[1].stop=30;
Method[1].maxit=2000;

Method[2].MethodName="Smart value function iteration";
Method[2].VI_PI=0;
Method[2].IP=0;
Method[2].BS=1;
Method[2].eps=0.01;
Method[2].stop=30;
Method[2].maxit=2000;


Method[3].MethodName="Full policy iteration";
Method[3].VI_PI=1;
Method[3].IP=0;
Method[3].BS=1;
Method[3].eps=0.01;
Method[3].stop=30;
Method[3].maxit=2000;


Method[4].MethodName="Partial policy function iteration";
Method[4].VI_PI=2;
Method[4].IP=0;
Method[4].VI_PI_K=30;
Method[4].BS=1;
Method[4].eps=0.01;
Method[4].stop=30;
Method[4].maxit=5000;

Method[5].MethodName="Smart value function iteration with interpolation";
Method[5].VI_PI=0;
Method[5].IP=1;
Method[5].BS=1;
Method[5].eps=0.000001;
Method[5].stop=30;
Method[5].maxit=10000;

Method[6].MethodName="Smart value function iteration with interpolation and partial policy function iteration";
Method[6].VI_PI=2;
Method[6].IP=2;
Method[6].VI_PI_K=5;
Method[6].BS=1;
Method[6].eps=0.01;
Method[6].stop=30;
Method[6].maxit=5000;

/* Parameters of the problem */
alpha=0.27;               @ elasticity of production with respect to capital                    @
beta=0.994;               @ discount factor                                                     @
delta=0.011;              @ rate of capital depreciation                                        @
eta=2.0;                  @ elasticity of marginal utility                                      @

/* Parameters of the algorithm */
kmin =0.75;               @ kmin*kstar is the lower bound of the grid @
kmax =1.25;               @ kmax*kstar is the upper bound of the grid @

/* select method number between 1 and 6 */
no=5;

_VI_eps=Method[no].eps;
_VI_stop=Method[no].stop;
_VI_BS=Method[no].BS;
_VI_Max=Method[no].maxit;
_VI_PI=Method[no].VI_PI;
_VI_PI_K=Method[no].VI_PI_K;
_IP=Method[no].IP;
_VI_BS=1;

/* Choose a vector of grid sizes */
@nvec={100,150,200,500,1000, 5000, 10000};@
nvec={250,500,1000,5000,10000};

clear _VI_grid, _VI_vec, _VI_der, _VI_kex; 

/* Parameters for the computation of Euler equation residuals */
kmin_e=0.8;               @ kmin_e*kstar is the lower bound @
kmax_e=1.2;               @ kmax_e*kstar is the upper bound @
nobs_e=200;               @ the number of residuals to be computed @

/* Name of output file */
output file=Ramsey2d4_4.txt reset;
MyDate;
?"Parameters of the model:";
?"alpha = " ftos(alpha,"*.*lf",7,3);
?"beta  = " ftos(beta, "*.*lf",7,3);
?"delta = " ftos(delta,"*.*lf",7,3);
?"eta   = " ftos(eta,  "*.*lf",7,3);
?"";
?"Algorithm: ";
?Method[no].Methodname;
?"";
?"Parameters of the algorithm:";
?"kmin  = " ftos(kmin, "*.*lf",7,3);
?"kmax  = " ftos(kmax, "*.*lf",7,3);
?"kmin_e= " ftos(kmin_e, "*.*lf",7,3);
?"kmax_e= " ftos(kmax_e, "*.*lf",7,3);
?"nobs_e= " ftos(nobs_e,"*.*lf",7,0);
?"";
?"eps   = " ftos(Method[no].eps,"*.*lf",10,8);
?"stop  = " ftos(Method[no].stop,"*.*lf",10,0);
?"IP    = " ftos(Method[no].IP,"*.*lf",10,0);
if Method[no].VI_PI==2;
?"k     = " ftos(Method[no].VI_PI_K,"*.*lf",10,0);
endif;
output off;

/* stationary solution */
kstar=((1- beta*(1-delta))/alpha*beta)^(1/(alpha-1));  @ stationary capital stock                          @
cstar=kstar^alpha - delta*kstar;                       @ stationary level of consumption                   @

kmin_g=kmin*kstar;
kmax_g=kmax*kstar;

kmin_e=kmin_e*kstar;
kmax_e=kmax_e*kstar;


evec=zeros(rows(nvec),1);
svec=evec;
lmax=rows(nvec);

policy=zeros(nobs_e,lmax);
eer=policy;
lmax=3;
for l (1,lmax,1);
    ngrid=nvec[l];

    /* Determination of grid */
    grid=seqa(kmin_g,(kmax_g-kmin_g)/(ngrid-1),ngrid);

    /* Initial value function */
    v0=rf(kstar,kstar)/(1-beta);
    v0=ones(ngrid,1).*v0;
    bsec=hsec;
    if no==1; x=Solve_VI_D1(&rf,beta,grid,v0); knext=grid[x]; endif;
    if no/=1;
        if _IP>0;
            knext=Solve_VI_D3(beta,grid,v0);
        else;
            x=Solve_VI_D2(&rf,beta,grid,v0);
            knext=grid[x];
        endif;
    endif;
    svec[l]=hsec-bsec;
    
    kvec=seqa(kmin_e,(kmax_e-kmin_e)/(nobs_e-1),nobs_e);
    
    dk=CSpline(grid,knext,1,0|0);
    
    for i (1,nobs_e,1);
        if _IP==2;
            k1=Splint(grid,knext,dk,1,kvec[i]);
            policy[i,l]=k1;
            k2=Splint(grid,knext,dk,1,k1);    
        else;
            k1=LIP(grid,knext,kvec[i]);
            policy[i,l]=k1;
            k2=LIP(grid,knext,k1);
        endif;
        c0=(kvec[i]^alpha)+(1-delta)*kvec[i] - k1;
        c1=(k1^alpha)+(1-delta)*k1-k2;
        rhs=beta*(c1^(-eta))*(1-delta+alpha*(k1^(alpha-1)));
        c1=rhs^(-1/eta);
        eer[i,l]=(c1/c0)-1; 
    endfor;
    
    evec[l]=maxc(abs(eer[.,l]));
    test="n= " $+ ftos(ngrid,"*.*lf",10,0) $+  " Run time= " $+ etstr(svec[l]) $+ " EER= " $+ ftos(evec[l],"*.*lf",10,8);
    output on;
    ?test;
    output off;
endfor;
save kvec, eer, policy;
end;

/* ------------------------------ Prodedures ------------------------------------ */

/* Definition of the utility function u=c^(1-eta)/(1-eta) as
** a function of this and next period's capital stock k1 and k2
** respectively */

proc(1)=rf(k1,k2);

  local c;
  c=k1^alpha + (1-delta)*k1 -k2;
  if c<0;
    c=miss(1,1);  @ set c to Gauss missing value code if c<0 @
  else;
    if eta==1;
      c=ln(c);
    else;  
      c=(c^(1-eta))/(1-eta);
    endif;
  endif;
 retp(c);
 
endp;
  

@  --------------------------  Solve_VI_D1 -------------------------------------

   Purpose: Computes the policy function for a simple deterministic, infinite
            horizon Ramsey model. In this model, capital K is the single
            variable) factor of production and depreciates at the constant
            rate delta (must be in the interval (0,1]). The production function
            f(K) is strictly increasing in K, the utility function u(C)
            is strictly concave in C=f(K)+(1-delta)*K-K'.
             
   Algorithm: The policy function is computed on a grid of n points,
              supplied by the user in the column vector grid. The
              method employed is iteration over the value function.


   Usage:  k=Solve_VI_D1(&rf,beta,grid,v0);

  Input :  &rf:= a pointer to the procedure that returns
                 u(f(K)+(1-delta)*K-K') for a pair (K,K') of points
                 from the grid.

          beta:= the discount factor from the agent's life-time utility function.

          grid:= n times 1 column vector with elements in ascending order that
                 represent valid stocks of the agent's capital stock.

          v0  := n times 1 column vector with the elements of the inital value function.
                  
  Output:  k  := n times 1 column vector. if i=1, ..., n is the
                 index of the a point from grid, k[i] is the index
                 of next period's optimal capital stock, i.e.,
                 if K'=g(K), g the policy function, then K'=grid[K[i]];

------------------------------------------------------------------------------------ @

  
proc(1)=Solve_VI_D1(&rf,beta,grid,v0);
                   
   local i, j, v1, v2, k, k1, k2, t, ngrid, dv, eps1, eps2, nc, di, w;
   local rf:proc;

/* Step 1: Initialize */
eps1=_VI_eps*(1-beta);      @ convergence criteria @
eps2=_VI_Stop;              @ stop after iterations with unchanged policy function @

nc=0;

ngrid=rows(grid);   @ number of grid points @
k=grid;           
k1=seqa(1,1,ngrid); @ intial policy function @
k2=zeros(ngrid,1);  
w=zeros(ngrid,1);

v1=v0;
v2=zeros(ngrid,1);
dv=1;

/* Step 2: Iterate over the value function */
DosWinOpen("Value Function Iteration",0|0|15|1|1);
cls;
t=1;
do until (t>_VI_Max) or (dv<eps1) or (nc>eps2);
  
  for i (1,ngrid,1);
  
    for j (1,ngrid,1);    
        w[j]=rf(k[i],k[j])+beta*v1[j];
    endfor;
    v2[i]=maxc(w);
    k2[i]=maxindc(w);
        
  endfor;
  di=sumc(k2 ./= k1);
  if di>=1; nc=0; else; nc=nc+1; endif;
  dv=maxc(maxc(abs(v2-v1))); 
  locate 5,5;
  ?"Iteration #= " ftos(t,"*.*lf",6,0);
  locate 6,5;
  ?"# of indices that have changed= " di;
  locate 7,5;
  ?"Largest element in v1-v0= " dv;
  locate 8,5;
  ?"# of consecutive iterations with constant policy function=" nc;
  v1=v2;
  k1=k2;
  t=t+1;
endo;
if t>_VI_Max;
   ?"Maximum number of iterations exceeded. Change _VI_Tmax!";
   ?"The computed solution may be inaccurate.";
   ?"Press any key...";wait;
endif;
if minc(k1)==1; locate 9,5; ?"Policy function hits lower bound of grid. Press any key to continue"; wait; endif;
if maxc(k1)==rows(k1); locate 10,5;?"Policy function hits upper bound of grid. Press any key to continue"; wait; endif;
retp(k1);

endp;


@  --------------------------  Solve_VI_D2 -------------------------------------

   Purpose: Computes the policy function for a simple deterministic, infinite
            horizon Ramsey model. In this model, capital K is the single
            variable) factor of production and depreciates at the constant
            rate delta (must be in the interval (0,1]). The production function
            f(K) is strictly increasing in K, the utility function u(C)
            is strictly concave in C=f(K)+(1-delta)*K-K'.
             
   Algorithm: The policy function is computed on a grid of n points,
              supplied by the user in the column vector grid. The
              method employed is iteration over the value function.

              Different from Solve_VI_D1 the algorithm makes use
              of the concavity of the value function and the montonicity
              of the policy function. Specifically, we use a binary
              search algorithm to locate the maximum on the rhs of
              the Bellman equation.


   Usage:  k=Solve_VI_D2(&rf,beta,grid,v0);

  Input :  &rf:= a pointer to the procedure that returns
                 u(f(K)+(1-delta)*K-K') for a pair (K,K') of points
                 from the grid.

          beta:= the discount factor from the agent's life-time utility function.

          grid:= n times 1 column vector with elements in ascending order that
                 represent valid stocks of the agent's capital stock.

          v0  := n times 1 column vector with the elements of the inital value function.
                  
  Output:  k  := n times 1 column vector. if i=1, ..., n is the
                 index of the a point from grid, k[i] is the index
                 of next period's optimal capital stock, i.e.,
                 if K'=g(K), g the policy function, then K'=grid[K[i]];

------------------------------------------------------------------------------------ @

  
proc(1)=Solve_VI_D2(&rf,beta,grid,v0);
                   
   local i, j, v1, v2, k, k1, k2, t, ngrid, dv, eps1, eps2, nc, di, w, js, jmin, jmax, jl, ju, amat, uvec;
   local rf:proc;

/* Step 1: Initialize */
eps1=_VI_eps*(1-beta);      @ convergence criteria @
eps2=_VI_Stop;              @ stop after iterations with unchanged policy function @

nc=0;

ngrid=rows(grid);   @ number of grid points @
k=grid;           
k1=seqa(1,1,ngrid); @ intial policy function @
k2=zeros(ngrid,1);  
w =zeros(3,1);

v1=v0;
v2=zeros(ngrid,1);
dv=1;

/* Step 2: Iterate over the value function */
DosWinOpen("Value Function Iteration",0|0|15|1|1);
cls;
t=1;
do until (t>_VI_Max) or (dv<eps1) or (nc>eps2);
  js=1;
  if _VI_PI==1; amat=zeros(ngrid,3); for j (1,ngrid,1); amat[j,1]=1; amat[j,2]=j; amat[j,3]=j; endfor; uvec=zeros(ngrid,1); endif; 
  if _VI_PI==2; amat=zeros(ngrid,3); uvec=zeros(ngrid,1); endif;

  for i (1,ngrid,1);  
    jmin=js;                
    jmax=ngrid;
    do while (jmax-jmin)>2;       @ the next lines implement the binary search algorithm @
        jl=floor((jmin+jmax)/2);
        ju=jl+1;
        w[1]=rf(k[i],k[jl])+beta*v1[jl];
        w[2]=rf(k[i],k[ju])+beta*v1[ju];
        if w[2]>w[1]; jmin=jl; else; jmax=ju; endif;
    endo;
    w[1]=rf(k[i],k[jmin])+beta*v1[jmin];
    if jmax>jmin;    w[2]=rf(k[i],k[jmin+1])+beta*v1[jmin+1]; else; w[2]=w[1]; endif;
    w[3]=rf(k[i],k[jmax])+beta*v1[jmax];
      js=maxindc(w);
   v2[i]=w[js]; 
      js=jmin+js-1;
   k2[i]=js;
   if _VI_PI==1; uvec[i]=rf(k[i],k[js]); if i==js; amat[i,1]=amat[i,1]-beta; else; amat=amat|(-beta~i~js); endif; endif;
   if _VI_PI==2; uvec[i]=rf(k[i],k[js]); amat[i,1]=1; amat[i,2]=i; amat[i,3]=js; endif;
   
  endfor;
  if _VI_PI==1; v2=SparseSolve(SparseFP(amat,ngrid,ngrid),uvec); endif;
  if _VI_PI==2; 
    for j (1,_VI_PI_K,1);           
          v2=uvec+beta*SparseTD(SparseFP(amat,ngrid,ngrid),v2);           
    endfor;
  endif;

  di=sumc(k2 ./= k1);
  if di>=1; nc=0; else; nc=nc+1; endif;
  dv=maxc(maxc(abs(v2-v1))); 
  locate 5,5;
  ?"Iteration #= " ftos(t,"*.*lf",6,0);
  locate 6,5;
  ?"# of indices that have changed= " di;
  locate 7,5;
  ?"Largest element in v1-v0= " dv;
  locate 8,5;
  ?"# of consecutive iterations with constant policy function=" nc;
  v1=v2;
  k1=k2;
  t=t+1;
endo;
if t>_VI_Max;
   ?"Maximum number of iterations exceeded. Change _VI_Tmax!";
   ?"The computed solution may be inaccurate.";
   ?"Press any key...";wait;
endif;
if minc(k1)==1; locate 9,5; ?"Policy function hits lower bound of grid"; endif;
if maxc(k1)==rows(k1); locate 10,5;?"Policy function hits upper bound of grid"; endif;
retp(k1);

endp;


@  --------------------------  Solve_VI_D3 -------------------------------------

   Purpose: Computes the policy function for a simple deterministic, infinite
            horizon Ramsey model. In this model, capital K is the single
            variable) factor of production and depreciates at the constant
            rate delta (must be in the interval (0,1]). The production function
            f(K) is strictly increasing in K, the utility function u(C)
            is strictly concave in C=f(K)+(1-delta)*K-K'.
             
   Algorithm: The policy function is computed on a grid of n points,
              supplied by the user in the column vector grid. The
              method employed is iteration over the value function.

              Different from Solve_VI_D1 the algorithm makes use
              of the concavity of the value function and the montonicity
              of the policy function. Specifically, we use a binary
              search algorithm to locate the maximum on the rhs of
              the Bellman equation.

              Different from Solve_VI_D2, the algorithm linearily 
              interpolates between grid points.


   Usage:  k=Solve_VI_D2(beta,grid,v0);

  Input : beta:= the discount factor from the agent's life-time utility function.

          grid:= n times 1 column vector with elements in ascending order that
                 represent valid stocks of the agent's capital stock.

          v0  := n times 1 column vector with the elements of the inital value function.
                  
  Output:  k  := n times 1 column vector.  if i=1, ..., n is the
                 index of the a point from grid, k[i] is the optimal
                 next period capital stock.


  Remarks: To be able to program this algorithm, the rhs of the bellman
           equation must be computed and the value function as a function
           of the next period state variable alone must be passed to
           an optimization routine. for this to work, the return function
           must be specified in rf(k0,k1) (no other name is allowed),
           and the global symbols _VI_vec, _VI_grid, and _VI_kex
           are needed.
           

------------------------------------------------------------------------------------ @

  
proc(1)=Solve_VI_D3(beta,grid,v0);
                   
   local i, j, v1, v2, k1, k2, t, ngrid, dv, eps1, di, w, js, jmin, jmax, jl, ju,
         ax, cx, bx;

   external matrix _VI_grid, _VI_vec, _VI_kex;

/* Step 1: Initialize */
eps1=_VI_eps*(1-beta);      @ convergence criteria @


ngrid=rows(grid);   @ number of grid points @

k1=grid;            @ intial policy function @
k2=zeros(ngrid,1);  
w =zeros(3,1);

v1=v0;
v2=zeros(ngrid,1);
dv=1;

_VI_vec=v0;
_VI_grid=grid;
if _IP==2; _VI_der=CSpline(grid,v0,1,0|0); endif;

/* Step 2: Iterate over the value function */
DosWinOpen("Value Function Iteration",0|0|15|1|1);
cls;
t=1;
do until (t>_VI_Max) or (dv<eps1);
  js=1;
  for i (1,ngrid,1);  
    _VI_kex=k1[i];  @ needed to compute the value function as a function of k alone @
    if _VI_BS;    
        jmin=js;                
        jmax=ngrid;
        do while (jmax-jmin)>2;       @ the next lines implement the binary search algorithm @
            jl=floor((jmin+jmax)/2);
            ju=jl+1;
            w[1]=rf(k1[i],k1[jl])+beta*v1[jl];
            w[2]=rf(k1[i],k1[ju])+beta*v1[ju];
            if w[2]>w[1]; jmin=jl; else; jmax=ju; endif;
        endo;
        w[1]=rf(k1[i],k1[jmin])+beta*v1[jmin];
        if jmax>jmin;    w[2]=rf(k1[i],k1[jmin+1])+beta*v1[jmin+1]; else; w[2]=w[1]; endif;
        w[3]=rf(k1[i],k1[jmax])+beta*v1[jmax];
        js=maxindc(w);
        js=jmin+js-1;
    else;    
        jmin=js;
        w[1]=rf(k1[i],k1[jmin])+beta*v1[jmin];         
        for jl (jmin+1,ngrid,1);                
             w[2]=rf(k1[i],k1[jl])+beta*v1[jl];
             if w[2]<=w[1]; js=jl-1; break; else; w[1]=w[2]; endif;                
        endfor; 
    endif;

    /* The next lines implement linear interpolation between grid points */
    if js==1;  /* boundary optimum, ax=bx=a[1]  */
                ax=k1[1];
                bx=ax+eps1*(k1[2]-k1[1]);
                cx=k1[2];
                if _rhs_bellman(k1[i],bx)<_rhs_bellman(k1[i],ax);
                    k2[i]=k1[1];
                else;                      
                    k2[i]=GSS(&_VI_valuefunction,k1[1],k1[2]);                    
                endif;
    elseif js==ngrid;   /* boundary optimum, bx=cx=a[n] */
                ax=k1[ngrid-1];
                cx=k1[ngrid];
                bx=cx-eps1*(k1[ngrid]-k1[ngrid-1]);
                if _rhs_bellman(k1[i],bx)<_rhs_bellman(k1[i],cx);
                    k2[i]=k1[ngrid];
                else;
                      k2[i]=GSS(&_VI_valuefunction,k1[ngrid-1],k1[ngrid]);                    
                endif;
    else;
                k2[i]=GSS(&_VI_valuefunction,k1[js-1],k1[js+1]);         
    endif;
    
    v2[i]=_rhs_bellman(k1[i],k2[i]);    

  endfor;
  if _VI_PI==2; for j (1,_VI_PI_K,1); _VI_vec=v2; for i (1,ngrid,1); v2[i]=_rhs_Bellman(k1[i],k2[i]); endfor; endfor; endif; 
  dv=maxc(maxc(abs(v2-v1))); 
  locate 5,5;
  ?"Iteration #= " ftos(t,"*.*lf",6,0);
  locate 6,5;
  ?"Largest element in v1-v0= " dv;
  
  _VI_vec=v2;
  if _IP==2; _VI_der=CSpline(grid,v2,1,0|0); endif;
  v1=v2;  
  t=t+1;
endo;
if t>_VI_Max;
   ?"Maximum number of iterations exceeded. Change _VI_Tmax!";
   ?"The computed solution may be inaccurate.";
   ?"Press any key...";wait;
endif;
if minc(k2)==k1[1];     locate 9,5; ?"Policy function hits lower bound of grid"; endif;
if maxc(k2)==k1[ngrid]; locate 10,5;?"Policy function hits upper bound of grid"; endif;

retp(k2);

endp;

/* vi=_rhs_bellman(x0,x1)
**
** returns the rhs of rf(x0,x1)+beta*v(x1) where v(x1) is found from
** linear interpolation between the to tabulated points in the vector
** _IV_vec that bracket k1. _VI_vec is an external global vector
** as well as _IV_grid which stores the grid point at which v(x) is
** tabulated
*/
proc(1)=_rhs_bellman(x0,x1);
  
 if _IP==1;
    retp(rf(x0,x1)+beta*LIP(_VI_grid,_VI_vec,x1));
 elseif _IP==2;
    retp(rf(x0,x1)+beta*Splint(_VI_grid,_VI_vec,_VI_der,1,x1));
 endif;  

endp;
 
/* vi=_VI_valuefunction(x1)
**
** returns the value function as a function of x1. Since the value function
** is also a function of the previous captial stock, the variable must
** be given in a gobal variable: _VI_kex.
*/

proc(1)=_VI_valuefunction(x1);

    retp(_rhs_bellman(_VI_kex,x1));

endp;

