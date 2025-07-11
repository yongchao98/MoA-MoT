import numpy as np
from scipy.optimize import minimize_scalar

def solve():
    """
    This function analyzes the logistic regression problem, computes the relevant
    parameters L and lambda, and determines the maximum learning rate for
    guaranteed convergence of gradient descent.
    """
    # Define a sample 1D non-separable dataset.
    # The points are non-separable because there's no w that can correctly
    # classify all points (e.g., if w>0, (-3,1) is misclassified; if w<0, (-2,-1) is misclassified).
    x = np.array([-3., -2., -1., 1., 2., 3.])
    y = np.array([1., -1., -1., 1., 1., -1.])
    N = len(x)

    print("Sample dataset (x, y):")
    for i in range(N):
        print(f"({x[i]:>4}, {y[i]:>2})")
    print("-" * 20)
    
    # Define the sigmoid function
    def sigma(t):
        t = np.clip(t, -500, 500) # for numerical stability
        return 1. / (1. + np.exp(-t))

    # Define the risk function R(w) and its second derivative R''(w)
    def R(w, x_data, y_data):
        z = y_data * x_data * w
        # Using the log-sum-exp trick for better numerical stability
        return np.mean(np.log(1 + np.exp(-z)))

    def R_double_prime(w, x_data, y_data):
        z = y_data * x_data * w
        s = sigma(z)
        return np.mean((y_data * x_data)**2 * s * (1 - s))

    # Find the optimal w* by minimizing the risk R(w)
    # For a 1D problem, minimize_scalar is efficient.
    res = minimize_scalar(lambda w: R(w, x, y), bounds=(-5, 5), method='bounded')
    w_star = res.x
    
    # Calculate L, the uniform smoothness constant
    # L = sup R''(w) = (1/4N) * sum(x_i^2)
    L = np.sum(x**2) / (4. * N)

    # Calculate lambda, the smoothness at the optimal point w*
    # lambda = R''(w*)
    lambda_val = R_double_prime(w_star, x, y)
    
    # The condition for guaranteed convergence is gamma < 2/L
    M = 2. / L

    print(f"The optimal solution is w* = {w_star:.4f}")
    print(f"The uniform smoothness constant is L = {L:.4f}")
    print(f"The smoothness at the optimum is lambda = {lambda_val:.4f}")
    print("\nBased on optimization theory for L-smooth convex functions,")
    print("the largest upper bound M for the learning rate gamma is 2/L.")
    print("\nThe final equation is:")
    print(f"M = 2 / L = 2 / {L:.4f} = {M:.4f}")

solve()
<<<C>>>