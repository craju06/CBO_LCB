print("MLCBBO Algorithm:")
import numpy as np
from scipy.stats import qmc, norm
from scipy.optimize import minimize
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
import warnings
import argparse
warnings.filterwarnings("ignore")
np.random.seed(123456)

#---------------------------------------------------------------------------
# Define the Objective and Constraint function:
def objective(arr):
    '''
    Define the objective function f(x).

    Parameters:
    - arr: 2-dimentional array.

    Returns:
    Return 1-dimentional function value.
    '''
    return -(np.cos((arr[0]-0.1)*arr[1]))**2- arr[0]*(np.sin(3*arr[0]+arr[1]))

def constraint(arr):
    '''
    Define the constraint function, c(x)<=0.

    Parameters:
    - arr: 2-dimentional array.

    Returns:
    Return 1-dimentional function value.
    '''
    t=np.arctan2(arr[0],(arr[1]+1E-9))
    return (arr[0]**2) + (arr[1]**2)-(2*np.cos(t)-0.5*np.cos(2*t)-0.25*np.cos(3*t)-0.125*np.cos(4*t))**2 -(2*np.sin(t))**2

#---------------------------------------------------------------------------
# Bounds
global l_bounds, u_bounds
l_bounds = [-2.25, -2.5] # lower bounds
u_bounds = [2.5, 1.75] # upper bounds

#---------------------------------------------------------------------------
# surrogate or approximation for the objective function
def surrogateF(model, X):
    """
    This function update the prior knowledge of objective function.

    Parameters:
    - model: Gaussian process (GP) surrogate model of objective function.
    - X: Sample prediction point.

    Returns:
    Returns the mean and standard deviation.
    """
    # catch any warning generated when making a prediction
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # return model prediction
        return model.predict(X, return_std=True)

#---------------------------------------------------------------------------
# surrogate or approximation for the constraint function
def surrogateC(model, Z):
    """
    This function update the prior knowledge of constraint function.

    Parameters:
    - model: Gaussian process (GP) surrogate model of constraint function.
    - X: Sample prediction point.

    Returns:
    Returns the mean and standard deviation.
    """
    # catch any warning generated when making a prediction
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # return model prediction
        return model.predict(Z, return_std=True)

#---------------------------------------------------------------------------
# Define MLCB acquistion function

def MLCB(Xsamples,f,modelf, modelc):
    '''
    Define the MLCB acquistion function.

    Parameters:
    - Xsamples: Sample points (2-dimensional array).
    - f: The True objective function value (1-dimensional array).
    - modelf: GP surrogate model of objective.
    - modelc: GP surrogate model of constraint.

    Returns:
    Return 1-dimensional array value of MLCB function.
    '''
    # prefict the mean and standard deviation at the sample point using the GP models.
    muf, stdf = surrogateF(modelf, Xsamples)
    muc, stdc = surrogateC(modelc, Xsamples)
    
    # Update the M value
    M = max(abs(f))
    # lambda parameter
    lambda1 = 1/norm.cdf(0.05)
    
    # Calculate the mean and standard deviation
    muF = (muf + M*(1-norm.cdf(-muc/(stdc+1E-9))))
    stdF = np.sqrt(stdf**2+(M**2)*(norm.cdf(-muc/(stdc+1E-9)))*(1-norm.cdf(-muc/(stdc+1E-9))))
    
    # Calculate the MLCB function
    Mlcb = muF - lambda1*stdF
    # Return the Mlcb value
    return Mlcb

#---------------------------------------------------------------------------
# Optimize MLCB acquistion function
def optacq(X,f,modelf,modelc):
    """
    Optimize the MLCB acquistion function

    Parameters:
    - Xsamples: Sample points (2-dimensional array).
    - f: The True objective function value (1-dimensional array).
    - modelf: GP surrogate model of objective.
    - modelc: GP surrogate model of constraint.

    Returns:
    Return the best sample point (2-dimensional array)
    """
    dim = np.shape(X)[1]
    
    # instant function to call the MLCB function
    def acq_fun(Xsamples):
        return MLCB(Xsamples.reshape(-1,dim),f,modelf,modelc)
    
    # Generate latin hypercube samples
    sampler = qmc.LatinHypercube(d=dim)
    sampleX = sampler.random(n=1000)
    Xsamples = qmc.scale(sampleX, l_bounds, u_bounds)
    
    # Evaluate the MLCB function at the samples
    scores = acq_fun(Xsamples)
    # Find the location of the minimum MLCB function value 
    ix = np.argmin(scores)
    
    # Store the minimum MLCB function value
    min_acq = scores.min()
    # Store the sample point of minimum function value in best_x
    best_x = Xsamples[ix]

    bounds = np.dstack((l_bounds,u_bounds))[0]
    
    # Run the L-BFGS-B algorithm as best_x is an initial input
    res = minimize(fun = acq_fun,
                   x0 = best_x,
                   bounds = bounds,
                   method = 'Nelder-Mead')
    if res.fun < min_acq:
        best_x = res.x
        min_acq = res.fun
    # Return the best sample point
    return best_x

#---------------------------------------------------------------------------
# Main MLCB-BO Algorithm

def MLCB_BO(X,f,c,K):
    """
    Run the M-LCBBO procedure

    Parameters:
    - X: Initial sample points (2-dimensional array).
    - f: Initial objective function values (1-dimensional array).
    - c: Initial constraint function values (1-dimensional array).
    - K: Number of additional function evaluation.

    Returns:
    Return the points (2-dimensional array) along with objective and constraint function values (1-dimensional array).
    """
    n_dims = np.shape(X)[1]
    kernel1 = Matern(length_scale=np.ones(n_dims),nu=2.5)
    # kernel1 = Matern(nu=2.5)
    modelf = GaussianProcessRegressor(kernel1)
    modelc = GaussianProcessRegressor(kernel1)
    modelf.fit(X,f)
    modelc.fit(X,c)
# main loop
    k = 0
    while (k<K):
        # print(k)
        x = optacq(X,f,modelf,modelc)
        f1 = objective(x)
        c1 = constraint(x)
        X = np.vstack((X,x))
        f = np.hstack((f,f1))
        c = np.hstack((c,c1))

        modelf.fit(X,f)
        modelc.fit(X,c)
        k= k+1
    return X,f,c

#---------------------------------------------------------------------------
# Simulation
# simu = 1 # Number of simulation
# n = 20   # Number of Initial sample 
# K = 100  # Number of addtional sample

parser = argparse.ArgumentParser()
parser.add_argument('arg1', type=int)
parser.add_argument('arg2', type=int)
parser.add_argument('arg3', type=int)
args = parser.parse_args()


n = args.arg1   # Number of Initial sample 
K = args.arg2  # Number of addtional sample
simu = args.arg3 # Number of simulation

XX =[]
ff = []
cc = []

for i in range(simu):
    print("Simulation:",i+1)
    sampler = qmc.LatinHypercube(d=2,seed = i+1)
    sampleX = sampler.random(n=n)
    X = qmc.scale(sampleX, l_bounds, u_bounds)
    f = np.asarray([objective(x) for x in X])
    c = np.asarray([constraint(z) for z in X])   
    x1,f1,c1 = MLCB_BO(X,f,c,K)
    XX.append(x1)
    ff.append(f1)
    cc.append(c1)

#---------------------------------------------------------------------------
# Calculate the best feasible objective value over the function evalutions
fx = np.copy(ff)
for i in range(np.shape(fx)[0]):
    idxx = np.where(cc[i]<=0)
    if idxx[0][0]!=0:
        fx[i][0:idxx[0][0]] = np.inf
    hx = np.zeros(np.shape(fx)[1])
    hx[idxx] = ff[i][idxx[0]]
    for j in range(np.shape(fx)[1]-1):
        if ff[i][j+1] < fx[i][j] and hx[j+1]!=0:
            fx[i][j+1] = ff[i][j+1]
        else:
            fx[i][j+1] = fx[i][j]

#---------------------------------------------------------------------------
# Calculate the average and quantiles
E2 = np.apply_along_axis(np.mean, 0, fx)
lower2 = np.apply_along_axis(np.quantile, 0, fx,0.05)
upper2 = np.apply_along_axis(np.quantile, 0, fx,0.95)

### print the results
# print(E2[45],E2[75],E2[119])
# print(lower2[45],lower2[75],lower2[119])
# print(upper2[45],upper2[75],upper2[119])

np.savetxt('MLCB.csv', E2, delimiter=',')

infeas_no = np.zeros(n+K)
for j in range(n+K):
    count = 0
    for i in range(simu):
        count += int(cc[i][j]>0)
    infeas_no[j] = count
np.savetxt('infeas_mlcb.csv', infeas_no, delimiter=',')

