import numpy as np
from scipy.optimize import root_scalar

def solve_and_explain():
    """
    This function demonstrates the convergence condition for gradient descent
    on a logistic regression problem.
    """
    # Step 1: Define a non-separable 1D dataset (x_i, y_i)
    # A dataset is non-separable if no w can satisfy sign(w*x_i) = y_i for all i.
    # For x_1=-2, y_1=1: sign(-2w) = 1 => w < 0
    # For x_2=-1, y_2=-1: sign(-w) = -1 => w > 0
    # These conditions conflict, so the dataset is non-separable, and a finite
    # optimal w_* exists.
    x = np.array([-2.0, -1.0, 1.0, 2.0])
    y = np.array([1.0, -1.0, 1.0, 1.0])
    N = len(x)

    print("--- Step 1: Problem Setup ---")
    print(f"Consider a non-separable dataset with N={N} points:")
    print(f"  x = {x}")
    print(f"  y = {y}")

    # Step 2: Define necessary functions based on the problem description
    def sigma(t):
        # The sigmoid function
        return 1.0 / (1.0 + np.exp(-t))

    def R_prime(w, x_arr, y_arr):
        # First derivative of the loss function R(w)
        # R'(w) = (1/N) * sum[ sigma(-y_i*w*x_i) * (-y_i*x_i) ]
        terms = sigma(-y_arr * w * x_arr) * (-y_arr * x_arr)
        return np.mean(terms)

    def R_double_prime(w, x_arr, y_arr):
        # Second derivative of the loss function R(w)
        # R''(w) = (1/N) * sum[ x_i^2 * sigma(y_i*w*x_i) * sigma(-y_i*w*x_i) ]
        t = y_arr * w * x_arr
        terms = (x_arr**2) * sigma(t) * sigma(-t)
        return np.mean(terms)

    # Step 3: Calculate L, the uniform (global) smoothness constant
    # L is the maximum value of R''(w), which occurs when the sigma product is 1/4.
    # L = (1/4N) * sum(x_i^2)
    sum_x_sq = np.sum(x**2)
    L = (1.0 / (4.0 * N)) * sum_x_sq

    print("\n--- Step 2: Calculate Global Smoothness L ---")
    print("The global smoothness constant L is the maximum curvature of the loss function.")
    print("The final equation for L is: L = (1 / (4 * N)) * sum(x_i^2)")
    print(f"  Using the numbers from our data:")
    print(f"  L = (1 / (4 * {N})) * {sum_x_sq} = {L:.4f}")

    # Step 4: Find the optimal w_* by solving R'(w_*) = 0
    # We use a numerical root finder on the first derivative.
    sol = root_scalar(lambda w: R_prime(w, x, y), bracket=[-10, 10], method='brentq')
    w_star = sol.root

    # Step 5: Calculate lambda, the smoothness (curvature) at the optimal point w_*
    lambda_val = R_double_prime(w_star, x, y)

    print("\n--- Step 3: Calculate Local Curvature lambda at the Optimum ---")
    print("First, we find the optimal weight w_* by solving R'(w_*) = 0.")
    print(f"  The optimal weight is w_* = {w_star:.4f}")
    print("Then, we calculate lambda = R''(w_*).")
    print(f"  lambda = {lambda_val:.4f}")
    print(f"As expected, L ({L:.4f}) >= lambda ({lambda_val:.4f}).")

    # Step 6: State the convergence condition and the final answer
    print("\n--- Step 4: Determine the Convergence Condition ---")
    print("For gradient descent to be guaranteed to converge from any starting point,")
    print("the learning rate gamma must be small enough to handle the steepest part of")
    print("the function, which is determined by the global smoothness L.")
    print("\nThe condition for guaranteed convergence is: gamma < 2 / L")
    
    bound_from_L = 2 / L
    bound_from_lambda = 2 / lambda_val

    print("\nThe final equation for the upper bound M is: M = 2 / L")
    print("Using the numbers from our calculation:")
    print(f"  M = 2 / {L:.4f} = {bound_from_L:.4f}")
    print("\nFor comparison, the bound calculated from lambda would be:")
    print(f"  2 / lambda = 2 / {lambda_val:.4f} = {bound_from_lambda:.4f}")
    print("Choosing a learning rate larger than 2/L (e.g., 3.0) could lead to divergence,")
    print("even if it is smaller than 2/lambda.")
    print("\nTherefore, the largest upper bound M that guarantees convergence for any")
    print("initialization and any data is 2/L.")

if __name__ == '__main__':
    solve_and_explain()
    print("\n<<<C>>>")
