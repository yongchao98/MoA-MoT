import numpy as np
from scipy.optimize import minimize_scalar

def solve():
    """
    This function demonstrates the calculation of the convergence bound M for gradient descent
    on a logistic regression problem.
    """
    # Step 1: Define a sample non-separable 1D dataset.
    # Let x_i = [1, 3, -2], y_i = [1, -1, -1].
    # The products y_i * x_i are [1, -3, 2].
    x = np.array([1., 3., -2.])
    y = np.array([1., -1., -1.])
    
    # Step 2: Define the loss function R(w) and its second derivative R''(w).
    # R(w) = (1/N) * sum(log(1 + exp(-y_i * w * x_i)))
    def R(w, x_data, y_data):
        z = y_data * w * x_data
        # Use log1p for numerical stability
        return np.mean(np.log1p(np.exp(-z)))

    # R''(w) = (1/N) * sum(x_i^2 * exp(-y_i*w*x_i) / (1 + exp(-y_i*w*x_i))^2)
    def R_double_prime(w, x_data, y_data):
        z = y_data * w * x_data
        exp_neg_z = np.exp(-z)
        return np.mean((x_data**2 * exp_neg_z) / (1 + exp_neg_z)**2)

    # Step 3: Calculate L, the uniform smoothness constant.
    # L is the supremum of R''(w). We find it by minimizing -R''(w).
    # The search space for the maximum of R''(w) is typically around w=0.
    res_L = minimize_scalar(lambda w: -R_double_prime(w, x, y), bounds=(-5, 5), method='bounded')
    L = -res_L.fun

    # Step 4: Calculate the optimal weight w* and the smoothness at the optimum, lambda.
    # Find w* by minimizing R(w).
    res_w_star = minimize_scalar(lambda w: R(w, x, y), bounds=(-5, 5), method='bounded')
    w_star = res_w_star.x
    
    # Lambda is the value of the second derivative at w*.
    lambda_val = R_double_prime(w_star, x, y)

    # Step 5: Determine the largest upper bound M.
    # The theory of convex optimization states that for guaranteed convergence, gamma < 2/L.
    # Thus, the largest upper bound M is 2/L.
    M_val = 2 / L

    # Print the results and the final equation.
    print("Based on optimization theory, the largest upper bound M for the learning rate gamma is given by the formula M = 2/L, where L is the uniform smoothness constant of the loss function.")
    print("\nFor the example dataset:")
    print(f"  x = {x}")
    print(f"  y = {y}")
    print("\nWe calculate the following values:")
    print(f"  L (supremum of R''(w)) = {L:.4f}")
    print(f"  w* (optimal weight) = {w_star:.4f}")
    print(f"  lambda (R''(w*)) = {lambda_val:.4f}")
    print("\nAs shown, L >= lambda. The convergence guarantee must use the global constant L.")
    print("\nThe final equation for M is:")
    print(f"  M = 2 / L")
    print(f"  M = 2 / {L:.4f}")
    print(f"  M = {M_val:.4f}")
    print("\nThis corresponds to answer choice C.")

solve()