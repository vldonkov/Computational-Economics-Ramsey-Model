The scripts are based on the book Dynamic General Equilibrium Modeling of Heer and Maussner, 2009.

The solution to a Stochastic Infinite-Horizon Ramsey Model via value function iterations
with/without linear interpolation on a discrete grid is contained in the script 'Main.m'.

The following scripts contain functions needed for running 'Main.m':
- 'SolveVIS.m':		computes the policy function for a Stochastic Infinite-Horizon Ramsey Model
- 'rf.m':		contains the one-period return function of the model
- 'MarkovAR.m':		approximates AR(1)-Process by Markov chain
- 'LIP.m':		linear interpolation function, implementing the interpolation between grid points
- 'VI_valuefunction.m': contains the value function as a function of x1 alone (linear interpolation)
- 'rhs_bellman.m':	returns the right-hand side of the Bellman equation (linear interpolation)
- 'GSS.m':		implementation of Golden Section Search, used to find the maximum of the interpolated value function
- 'BLIP.m':		bilinear intepolation, used for the computation of the policy function
- 'Euler.m':		calculates the Euler residuals
- 'GetRhs.m':		right-hand side of Euler equation, used in 'Euler.m'
- 'PF.m':		calculates the policy function via bilinear interpolation, used for the calculation of the Euler residuals
- 'GH_INT4.m':		Gauss-Hermite integration over one-dimensional space, used in 'Euler.m'
- 'MachEps.m':		computes the machine epsilon