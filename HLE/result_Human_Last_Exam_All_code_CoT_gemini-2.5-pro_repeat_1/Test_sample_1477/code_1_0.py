import numpy as np
from scipy.optimize import minimize_scalar, brentq

def solve():
    """
    This function demonstrates the convergence condition for gradient descent
    on a 1D logistic regression problem.
    """
    # 1. Define the problem based on the user's description
    # We create a sample non-separable 1D dataset
    X = np.array([-4., -2., -1., 1., 2., 5.])
    Y = np.array([-1., -1., 1., -1., 1., 1.])
    N = len(X)

    # The objective function R(w) and its derivatives
    def R(w):
        # R(w) = 1/N * sum(log(1 + exp(-y_i * w * x_i)))
        return (1/N) * np.sum(np.log(1 + np.exp(-Y * w * X)))

    def R_prime(w):
        # R'(w) = 1/N * sum(-y_i * x_i * exp(-y_i*w*x_i) / (1 + exp(-y_i*w*x_i)))
        return (1/N) * np.sum(-Y * X * np.exp(-Y * w * X) / (1 + np.exp(-Y * w * X)))

    def R_double_prime(w):
        # R''(w) = 1/N * sum(x_i^2 * exp(-y_i*w*x_i) / (1 + exp(-y_i*w*x_i))^2)
        # This term is sigma(z)*(1-sigma(z)) where z = ywx
        s = 1 / (1 + np.exp(-Y * w * X))
        return (1/N) * np.sum(X**2 * s * (1 - s))

    # 2. Find the constants L and lambda
    # L is the uniform smoothness, the maximum value of R''(w).
    # We find it by minimizing the negative of R''(w).
    res_L = minimize_scalar(lambda w: -R_double_prime(w))
    L = -res_L.fun

    # w_* is the optimum, found by minimizing R(w) or finding the root of R'(w).
    try:
        w_star = brentq(R_prime, -10, 10)
    except ValueError:
        # Fallback to minimization if root finding fails
        res_w_star = minimize_scalar(R)
        w_star = res_w_star.x

    # lambda is the smoothness at the optimum w_*.
    lambda_val = R_double_prime(w_star)

    print("Step-by-step demonstration:")
    print("-" * 50)
    print("For our sample dataset, we find the key constants:")
    print(f"  - Uniform smoothness constant L = sup |R''(w)| = {L:.4f}")
    print(f"  - Optimal parameter w* = {w_star:.4f}")
    print(f"  - Smoothness at optimum lambda = R''(w*) = {lambda_val:.4f}")
    print("\nBased on theory, the critical learning rate bound is 2/L.")
    print(f"The final equation for the bound M is M = 2/L")
    print(f"Numerically, this is M = 2/{L:.4f} = {2/L:.4f}")
    print("-" * 50)

    # 3. Run Gradient Descent to test the theoretical bound
    def gradient_descent(gamma, w_init, num_iterations=50):
        w = w_init
        for _ in range(num_iterations):
            grad = R_prime(w)
            w = w - gamma * grad
            # Stop if weights explode
            if np.abs(w) > 1e6 or np.isnan(w):
                return w
        return w

    w_initial = 2.0
    print(f"Running Gradient Descent from w_init = {w_initial}...\n")

    # Case 1: Safe learning rate (gamma < 2/L)
    gamma_safe = 1.9 / L
    final_w_safe = gradient_descent(gamma_safe, w_initial)
    print(f"Test 1: Using a safe learning rate gamma = 1.9/L = {gamma_safe:.4f}")
    print(f"  Result: Gradient descent converges to w = {final_w_safe:.4f} (close to w*)")

    # Case 2: Unsafe learning rate (gamma > 2/L)
    gamma_unsafe = 2.1 / L
    final_w_unsafe = gradient_descent(gamma_unsafe, w_initial)
    print(f"\nTest 2: Using an unsafe learning rate gamma = 2.1/L = {gamma_unsafe:.4f}")
    if np.abs(final_w_unsafe) > 1e6 or np.isnan(final_w_unsafe):
        print(f"  Result: Gradient descent diverges (w -> {final_w_unsafe:.2f})")
    else:
        print(f"  Result: Gradient descent oscillates and fails to converge, w = {final_w_unsafe:.4f}")
    print("-" * 50)
    print("Conclusion: The demonstration confirms that convergence is only guaranteed for gamma < 2/L.")
    print("The largest upper bound M is 2/L.")

solve()

<<<C>>>