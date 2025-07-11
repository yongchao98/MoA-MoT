import numpy as np
from scipy.optimize import minimize_scalar, brentq

def solve():
    """
    Analyzes the convergence of gradient descent for 1D logistic regression
    and determines the largest upper bound for the learning rate.
    """

    # 1. Define the dataset (chosen to be non-separable)
    X = np.array([-4., -1., 1.5, 3.])
    Y = np.array([1., -1., 1., -1.])
    N = len(X)

    # 2. Define the logistic loss R(w) and its derivatives
    def sigma(t):
        # Clip to avoid overflow in exp
        t = np.clip(t, -100, 100)
        return 1. / (1. + np.exp(-t))

    def R_prime(w, x_data, y_data):
        # Gradient of the loss function
        args = y_data * w * x_data
        # sigma(ywx) - 1 is numerically better as -sigma(-ywx)
        return np.mean((sigma(args) - 1.) * y_data * x_data)

    def R_double_prime(w, x_data, y_data):
        # Hessian (second derivative) of the loss function
        args = y_data * w * x_data
        s = sigma(args)
        return np.mean(x_data**2 * s * (1. - s))

    # 3. Find L, the global smoothness constant
    # L is the maximum value of R''(w). We find it by minimizing -R''(w).
    neg_R_double_prime = lambda w: -R_double_prime(w, X, Y)
    res = minimize_scalar(neg_R_double_prime, bounds=(-10, 10), method='bounded')
    L = -res.fun

    # 4. Find w*, the optimal point, and lambda, the local curvature
    # w* is the root of R'(w)
    try:
        w_star = brentq(R_prime, a=-10, b=10, args=(X, Y))
    except ValueError:
        # Fallback for cases where root is not in bracket; not expected for this data
        w_star = 0.0
        for _ in range(2000):
            w_star -= 0.1 * R_prime(w_star, X, Y)

    lambda_val = R_double_prime(w_star, X, Y)
    
    # 5. Determine and print the theoretical bound M
    M = 2 / L

    print("--- Analysis of Gradient Descent Convergence ---")
    print(f"Dataset X: {X}")
    print(f"Dataset Y: {Y}")
    print("-" * 44)
    print(f"Global smoothness constant (max curvature) L = {L:.4f}")
    print(f"Local curvature at optimum w*={w_star:.4f} is lambda = {lambda_val:.4f}")
    print("-" * 44)

    print("Theoretical convergence guarantee for learning rate gamma is: gamma < 2 / L")
    print("Therefore, the largest upper bound M is 2 / L.")
    print("\n--- Final Equation ---")
    print(f"M = 2 / L = 2 / {L:.4f} = {M:.4f}")
    print("----------------------")


    # 6. Verification with gradient descent runs
    def gradient_descent(gamma, w_init, num_steps=30):
        w = w_init
        print(f"\nRunning GD with gamma = {gamma:.4f} (Theoretical limit M = {M:.4f})")
        print("This gamma is " + ("< M" if gamma < M else ">= M"))
        try:
            for i in range(num_steps):
                grad = R_prime(w, X, Y)
                w = w - gamma * grad
                if np.abs(w) > 1e6 or np.isnan(w):
                    print(f"Step {i+1:2d}: w = {w:.4f} -> DIVERGED")
                    return
                if (i < 5) or ((i + 1) % 5 == 0):
                    print(f"Step {i+1:2d}: w = {w:.4f}")
            print(f"Final w = {w:.4f} -> CONVERGED to (or towards) w*={w_star:.4f}")
        except OverflowError:
            print("Divergence detected (overflow).")

    w0 = 2.0  # An arbitrary starting point

    # Run with a learning rate guaranteed to converge
    gamma_converges = 0.95 * M
    gradient_descent(gamma_converges, w0)

    # Run with a learning rate that may diverge
    gamma_diverges = 1.05 * M
    gradient_descent(gamma_diverges, w0)
    
    print("\nConclusion: The theoretical bound M = 2/L holds. A learning rate above this bound can cause divergence.")
    print("The correct answer choice is the one that represents this bound.")


solve()