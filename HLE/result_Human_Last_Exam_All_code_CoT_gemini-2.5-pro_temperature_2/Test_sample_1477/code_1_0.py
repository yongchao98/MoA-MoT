import numpy as np
from scipy.optimize import root_scalar

def solve_logistic_regression_bound():
    """
    Analyzes the convergence bound for 1D logistic regression.

    This function follows the plan:
    1. Define a sample non-separable 1D dataset.
    2. Define the necessary functions: sigmoid, risk derivative R'(w), and second derivative R''(w).
    3. Calculate the uniform smoothness L.
    4. Numerically find the optimal weight w_*.
    5. Calculate the local smoothness at the optimum, lambda = R''(w_*).
    6. Compute the values corresponding to the answer choices and print the result.
    """

    # 1. Define a sample non-separable dataset (N=2)
    # "Non-separable" in 1D means that the signs of y_i * x_i are not all the same.
    x = np.array([2.0, 1.0])
    y = np.array([1.0, -1.0])
    N = len(x)

    # 2. Define helper functions
    def sigmoid(t):
        return 1 / (1 + np.exp(-t))

    def risk_prime(w, x_data, y_data):
        # R'(w) = -1/N * sum(y_i * x_i * sigma(-y_i * w * x_i))
        t = -y_data * w * x_data
        return -np.mean(y_data * x_data * sigmoid(t))

    def risk_double_prime(w, x_data, y_data):
        # R''(w) = 1/N * sum(x_i^2 * sigma(-y_i*w*x_i) * (1-sigma(-y_i*w*x_i)))
        t = -y_data * w * x_data
        sig_t = sigmoid(t)
        return np.mean(x_data**2 * sig_t * (1 - sig_t))

    # 3. Calculate L
    # L = sup_w R''(w). The max of sigma(t)*(1-sigma(t)) is 1/4 at t=0.
    # So L is R''(0).
    L = np.mean(x**2 * 0.25)
    
    # 4. Find the optimal weight w_*, where R'(w_*) = 0
    # We use a numerical root finder on the R'(w) function.
    sol = root_scalar(lambda w: risk_prime(w, x, y), bracket=[-5, 5], method='brentq')
    w_star = sol.root

    # 5. Calculate lambda = R''(w_*)
    lambda_val = risk_double_prime(w_star, x, y)
    
    # 6. Compute and display the results for the answer choices
    print("Based on the provided data:")
    print(f"  x = {x}, y = {y}")
    print("-" * 30)
    print("Step 1: Calculate the uniform smoothness constant L.")
    # The term sigma(t)*(1-sigma(t)) has a maximum value of 0.25.
    # L = sup_w R''(w) = (1/N) * sum(x_i^2 * max(sigma(t)*(1-sigma(t))))
    L_calc_str = f"L = (1/{N}) * ({x[0]**2} + {x[1]**2}) * 0.25 = {L}"
    print(L_calc_str)
    print("-" * 30)

    print(f"Step 2: Find the optimal weight w_* where R'(w_*) = 0.")
    print(f"  Numerically found w_* = {w_star:.4f}")
    print("\nStep 3: Calculate lambda = R''(w_*).")
    print(f"  lambda = R''({w_star:.4f}) = {lambda_val:.4f}")
    print("-" * 30)
    
    print("Step 4: Evaluate potential bounds M based on L and lambda.")
    
    # Option A
    m_a = 1 / L
    print(f"A. 1/L = 1/{L:.4f} = {m_a:.4f}")
    
    # Option B
    m_b = 2 / lambda_val
    print(f"B. 2/lambda = 2/{lambda_val:.4f} = {m_b:.4f}")
    
    # Option C
    m_c = 2 / L
    print(f"C. 2/L = 2/{L:.4f} = {m_c:.4f}")

    # Option D
    print(f"D. infinity")
    
    # Option E
    m_e = max(1/lambda_val, 2/L)
    print(f"E. max(1/lambda, 2/L) = max({1/lambda_val:.4f}, {2/L:.4f}) = {m_e:.4f}")
    
    print("-" * 30)
    print("The theoretical condition for guaranteed global convergence of gradient descent")
    print("for an L-smooth convex function is that the learning rate gamma must be < 2/L.")
    print("Therefore, the largest upper bound M is 2/L.")

solve_logistic_regression_bound()
<<<C>>>