import numpy as np
from scipy.optimize import fsolve, minimize_scalar

def solve_and_demonstrate():
    """
    This function sets up a 1D logistic regression problem, calculates L and lambda,
    and demonstrates the convergence condition for gradient descent.
    """
    # Step 1: Define a non-separable 1D dataset
    X = np.array([1.0, 2.0])
    Y = np.array([1.0, -1.0])
    N = len(X)
    print(f"Using non-separable dataset (x, y): {list(zip(X, Y))}\n")

    # Step 2: Define the loss function R(w) and its derivatives
    def R(w, X, Y):
        z = Y * w * X
        return np.mean(np.log(1 + np.exp(-z)))

    def R_grad(w, X, Y):
        z = Y * w * X
        sigma_z = 1 / (1 + np.exp(-z))
        return np.mean(Y * X * (sigma_z - 1))

    def R_hess(w, X, Y):
        z = Y * w * X
        # equivalent to x^2 * sigma(z) * sigma(-z)
        return np.mean((X**2) * np.exp(-z) / (1 + np.exp(-z))**2)

    # Step 3: Find the optimal weight w_star by finding the root of the gradient
    # The fsolve function finds the roots of a function.
    w_star_solution = fsolve(lambda w: R_grad(w, X, Y), x0=0.0)
    w_star = w_star_solution[0]
    print(f"Optimal weight w_* found at: {w_star:.4f}")

    # Step 4: Find the smoothness constant L = sup_w R''(w)
    # We maximize R_hess by minimizing its negative.
    res = minimize_scalar(lambda w: -R_hess(w, X, Y))
    w_max_hess = res.x
    L = -res.fun
    print(f"Smoothness constant L (max of R''(w)) = {L:.4f}, found at w = {w_max_hess:.4f}")

    # Step 5: Calculate lambda = R''(w_*)
    lambda_val = R_hess(w_star, X, Y)
    print(f"Curvature at optimum lambda (R''(w*)) = {lambda_val:.4f}\n")

    # The theoretical result for the largest upper bound M
    M = 2 / L
    
    print("--- Theoretical Analysis ---")
    print("The convergence of gradient descent is guaranteed for a learning rate gamma < 2/L.")
    print("The largest upper bound M for gamma is therefore 2/L.")
    print("\n--- Final Equation ---")
    print(f"The equation for the bound is M = 2 / L")
    print(f"Using the numbers from our example:")
    print(f"numerator = 2")
    print(f"denominator = L = {L:.4f}")
    print(f"M = 2 / {L:.4f} = {M:.4f}")
    print("------------------------\n")


    # Step 6: Demonstrate Gradient Descent
    def gradient_descent(gamma, w_init=2.0, num_iterations=50):
        w = w_init
        print(f"--- Running GD with gamma = {gamma:.4f} ---")
        print(f"Starting at w = {w_init:.4f}")
        for i in range(num_iterations):
            grad = R_grad(w, X, Y)
            w = w - gamma * grad
            if i < 5 or i > num_iterations - 6 :
                print(f"Iter {i+1:2d}: w = {w:.6f}")
            if np.abs(w) > 1e6 or np.isnan(w):
              print("!!! Divergence detected !!!")
              return
        print(f"Final w after {num_iterations} iterations: {w:.6f}")
        print(f"Converged towards w_* = {w_star:.6f}")


    # Demonstrate convergence for gamma < 2/L
    gamma_good = 0.95 * M
    gradient_descent(gamma_good)

    print("\n")

    # Demonstrate divergence for gamma > 2/L
    gamma_bad = 1.05 * M
    gradient_descent(gamma_bad)

solve_and_demonstrate()