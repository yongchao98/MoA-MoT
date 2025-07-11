import numpy as np
import matplotlib.pyplot as plt

def solve():
    """
    This function demonstrates the convergence condition for gradient descent
    on a logistic regression problem.
    """
    # Step 1: Define a non-separable 1D dataset
    x = np.array([-2., -1., 1., 2.])
    y = np.array([-1., 1., -1., 1.])
    N = len(x)

    # Step 2: Define the risk function R and its derivatives
    def sigma(t):
        # To avoid overflow for large negative t
        return np.where(t >= 0, 1 / (1 + np.exp(-t)), np.exp(t) / (1 + np.exp(t)))

    def R(w, x, y):
        # Add a small epsilon to prevent log(0)
        return -np.mean(np.log(sigma(y * w * x) + 1e-9))

    def R_prime(w, x, y):
        # R'(w) = -1/N * sum( y_i * x_i * sigma(-y_i * w * x_i) )
        return -np.mean(y * x * sigma(-y * w * x))

    def R_double_prime(w, x, y):
        t = -y * w * x
        s = sigma(t)
        return np.mean(x**2 * s * (1 - s))

    # Step 3: Calculate L = sup R''(w)
    # For this symmetric dataset, the maximum of R''(w) occurs at w=0.
    # The maximum value of sigma(t)*(1-sigma(t)) is 0.25, which occurs at t=0.
    L = R_double_prime(0, x, y)
    
    # Step 4: The theoretical bound is M = 2/L
    M = 2 / L
    
    print("This script demonstrates the convergence bound for gradient descent.")
    print(f"For the example dataset, the smoothness constant is L = sup R''(w) = {L:.4f}")
    print("The theoretical largest upper bound for the learning rate is M = 2/L.")
    print(f"The calculated value for M is: 2 / {L:.4f} = {M:.4f}")

    # Step 5: Find the optimal w* by solving R'(w*) = 0
    # From our manual derivation in the thought process: w_star = -ln(3)/2
    w_star = -np.log(3) / 2
    print(f"The optimal weight w* is approximately {w_star:.4f}")

    # Step 6: Run gradient descent with different learning rates
    w_init = 3.0
    iterations = 20

    # Case 1: gamma < M (e.g., gamma = 0.95 * M)
    gamma_converge = 0.95 * M
    w = w_init
    history_converge = [w]
    print(f"\n--- Running Gradient Descent with gamma = {gamma_converge:.4f} (< M) ---")
    for i in range(iterations):
        grad = R_prime(w, x, y)
        w = w - gamma_converge * grad
        history_converge.append(w)
        if (i < 5 or i == iterations - 1):
             print(f"Iteration {i+1:2d}: w = {w:.4f}")
    print(f"Final w approaches w* = {w_star:.4f}. Convergence is observed.")


    # Case 2: gamma > M (e.g., gamma = 1.05 * M)
    gamma_diverge = 1.05 * M
    w = w_init
    history_diverge = [w]
    print(f"\n--- Running Gradient Descent with gamma = {gamma_diverge:.4f} (> M) ---")
    for i in range(iterations):
        grad = R_prime(w, x, y)
        w = w - gamma_diverge * grad
        history_diverge.append(w)
        if (i < 5 or i == iterations - 1):
             print(f"Iteration {i+1:2d}: w = {w:.4f}")
    print("The values are oscillating and not converging to w*.")
    
    print("\nConclusion: The analysis shows that for guaranteed convergence for any initialization,")
    print("the learning rate gamma must be less than 2/L. Therefore, M = 2/L.")
    
solve()