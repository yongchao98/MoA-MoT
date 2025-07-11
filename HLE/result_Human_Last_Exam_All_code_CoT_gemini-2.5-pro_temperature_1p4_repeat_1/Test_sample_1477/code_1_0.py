import numpy as np
from scipy.optimize import minimize
import warnings

def solve():
    """
    This script solves for the largest upper bound M for the learning rate
    in logistic regression gradient descent and demonstrates the result.
    """

    # 1. Define a non-separable dataset (x_i, y_i)
    # The points are not linearly separable by a plane through the origin.
    x = np.array([-3.0, -1.0, 1.0, 2.0])
    y = np.array([1.0, -1.0, -1.0, 1.0])
    N = len(x)

    # 2. Define sigma, the loss function R(w), and its derivatives
    def sigma(t):
        # Clip to avoid overflow in exp, which can cause RuntimeWarning
        t = np.clip(t, -500, 500)
        return 1.0 / (1.0 + np.exp(-t))

    def R(w, x_data, y_data):
        # The loss function. w is a single value, but scipy expects an array.
        w_val = w[0] if isinstance(w, (np.ndarray, list)) else w
        # Clip sigma's output to avoid log(0)
        sig_vals = sigma(y_data * w_val * x_data)
        sig_vals = np.clip(sig_vals, 1e-15, 1 - 1e-15)
        return -np.mean(np.log(sig_vals))

    def R_prime(w, x_data, y_data):
        # First derivative of R(w)
        return -np.mean(y_data * x_data * sigma(-y_data * w * x_data))

    def R_double_prime(w, x_data, y_data):
        # Second derivative of R(w)
        u = y_data * w * x_data
        sig_u = sigma(u)
        return np.mean(x_data**2 * sig_u * (1 - sig_u))

    # 3. Calculate L, the uniform smoothness constant
    # L is the maximum value of R''(w). The term sigma(t)*(1-sigma(t)) is
    # maximized at t=0, with a value of 0.25.
    L = (1 / (4.0 * N)) * np.sum(x**2)
    
    # 4. Find the optimal w_star by minimizing R(w)
    result = minimize(lambda w: R(w, x, y), x0=np.array([0.0]), method='BFGS')
    w_star = result.x[0]

    # 5. Calculate lambda, the smoothness at the optimum w_star
    lambda_val = R_double_prime(w_star, x, y)

    # 6. Calculate the largest upper bound M for the learning rate
    M = 2.0 / L

    print("--- Problem Setup ---")
    print(f"Data points (x, y): {list(zip(x, y))}")
    print(f"Number of points N = {N}")
    print(f"Sum of x_i^2 = {np.sum(x**2):.4f}\n")

    print("--- Smoothness Calculation ---")
    print("L is the uniform smoothness constant (max curvature).")
    print(f"L = (1 / (4 * N)) * sum(x_i^2) = (1 / (4 * {N})) * {np.sum(x**2):.2f} = {L:.4f}")
    print("lambda is the smoothness at the optimal point w* (local curvature).")
    print(f"Found optimal w* = {w_star:.4f}")
    print(f"lambda = R''(w*) = {lambda_val:.4f}\n")
    
    print("--- Convergence Bound Calculation ---")
    print("The convergence of gradient descent is guaranteed for a learning rate gamma < 2/L.")
    print("This is because the guarantee must hold for any starting point, including")
    print("regions of maximum curvature L.")
    print(f"The largest upper bound is M = 2 / L")
    print(f"M = 2 / {L:.4f} = {M:.4f}\n")
    
    # 7. Demonstrate convergence vs. divergence
    print("--- Demonstration ---")
    def gradient_descent(gamma, w_init, iterations):
        w = w_init
        for _ in range(iterations):
            grad = R_prime(w, x, y)
            w = w - gamma * grad
        return w

    w_init_demo = 5.0
    iterations = 50

    # Case 1: gamma < M (should converge)
    gamma_converge = M * 0.99
    w_final_converge = gradient_descent(gamma_converge, w_init_demo, iterations)
    print(f"Running GD with gamma = {gamma_converge:.4f} (< M)...")
    print(f"w starts at {w_init_demo} and converges to {w_final_converge:.4f} (w* is {w_star:.4f}).\n")

    # Case 2: gamma > M (should diverge)
    # Catch warnings from overflow in exp
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        try:
            gamma_diverge = M * 1.01
            w_final_diverge = gradient_descent(gamma_diverge, w_init_demo, iterations)
            print(f"Running GD with gamma = {gamma_diverge:.4f} (> M)...")
            if np.isnan(w_final_diverge) or np.isinf(w_final_diverge):
                print(f"w starts at {w_init_demo} and diverges to {w_final_diverge}.")
            else:
                 print(f"w starts at {w_init_demo} and moves to {w_final_diverge} (diverging behavior).")
        except RuntimeWarning:
            print(f"Running GD with gamma = {gamma_diverge:.4f} (> M)...")
            print(f"w starts at {w_init_demo} and diverges (OverflowError encountered).")
            
solve()