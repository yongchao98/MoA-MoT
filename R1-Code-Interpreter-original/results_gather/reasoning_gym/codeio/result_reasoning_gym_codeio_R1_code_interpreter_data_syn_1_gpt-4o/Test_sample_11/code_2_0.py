import numpy as np
from numpy import linalg as la
import math

# Hard thresholding function
def hardThreshold(x, k):
    p = x.shape[0]
    t = np.sort(np.abs(x))[::-1]    
    threshold = t[k-1]
    j = (np.abs(x) < threshold)
    x[j] = 0
    return x

# IHT algorithm
def IHT(y, A, k, iters, epsilon):
    p = A.shape[1]   # Length of original signal    
    n = A.shape[0]   # Length of measurement vector
    
    x_new = np.zeros(p)    # Initial estimate    
    At = np.transpose(A)   # Transpose of A

    for i in range(iters):
        x_old = x_new
    
        # Compute gradient
        grad = -At @ (y - A @ x_new)
    
        # Perform gradient step
        x_temp = x_old - 0.5 * grad    
    
        # Perform hard thresholding step
        x_new = hardThreshold(x_temp, k)
    
        if (la.norm(x_new - x_old, 2) / la.norm(x_new, 2)) < epsilon:
            break
    
    return x_new

# main function
def main_solution(n, p, k, iters, epsilon):
    # Generate a p-dimensional zero vector
    x_star = np.zeros(p)
    # Randomly sample k indices in the range [1:p]
    x_star_ind = np.random.choice(p, k, replace=False) 
    # Set x_star_ind with k random elements from Gaussian distribution
    x_star[x_star_ind] = np.random.randn(k)
    # Normalize
    x_star = (1 / la.norm(x_star, 2)) * x_star

    # Generate sensing matrix
    A = (1 / math.sqrt(n)) * np.random.randn(n, p)

    # Observation model
    y = A @ x_star

    # Run IHT algorithm
    x_IHT = IHT(y, A, k, iters, epsilon)

    # Return the norm of the difference between the true and estimated vectors
    return la.norm(x_star - x_IHT, 2)

# Test different parameters to find a match
n = 80
p = 200
k = 20
iters = 500
epsilon = 1e-4

result = main_solution(n, p, k, iters, epsilon)
print(result)