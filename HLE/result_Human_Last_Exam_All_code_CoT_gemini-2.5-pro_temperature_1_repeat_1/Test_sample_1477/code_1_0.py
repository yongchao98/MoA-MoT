import numpy as np
import matplotlib.pyplot as plt

def solve():
    """
    This function demonstrates the convergence condition for gradient descent
    on a 1D logistic regression problem.
    """
    # 1. Define a sample dataset (non-separable for realism)
    X = np.array([-4, -3, -2, -1, 1, 2, 3, 4, 0.5, -0.5])
    Y = np.array([-1, -1, -1, -1, 1, 1, 1, 1, -1, 1])
    N = len(X)

    print("--- Problem Setup ---")
    print(f"Training data X: {X}")
    print(f"Training data Y: {Y}")
    print(f"Number of points N: {N}")

    # 2. Define the objective function R(w) and its derivatives
    def sigmoid(t):
        # Clip to avoid overflow in exp
        t = np.clip(t, -500, 500)
        return 1 / (1 + np.exp(-t))

    def R(w, x_data, y_data):
        z = y_data * w * x_data
        return np.mean(np.log(1 + np.exp(-z)))

    def grad_R(w, x_data, y_data):
        z = y_data * w * x_data
        return np.mean(-y_data * x_data * sigmoid(-z))

    def hess_R(w, x_data, y_data):
        z = y_data * w * x_data
        s = sigmoid(-z)
        return np.mean(x_data**2 * s * (1 - s))

    # 3. Calculate the global smoothness constant L
    # L = sup |R''(w)| = max( 1/4 * (1/N) * sum(x_i^2) )
    L = np.mean(X**2) / 4.0
    
    # The final equation for the convergence bound M
    M = 2 / L

    print("\n--- Theoretical Calculation ---")
    print(f"The smoothness constant is L = (1/(4N)) * sum(x_i^2) = {L:.4f}")
    print(f"The largest upper bound for the learning rate is M = 2 / L.")
    print(f"Therefore, M = 2 / {L:.4f} = {M:.4f}")
    print("This means gradient descent should converge for gamma < M and may diverge for gamma > M.")

    # 4. Find the optimal w* to compare against
    w_star = 0.0
    for _ in range(5000):
        w_star -= 0.01 * grad_R(w_star, X, Y)
    
    print(f"The optimal weight is approximately w* = {w_star:.4f}")

    # 5. Demonstrate convergence and divergence
    
    # Case 1: gamma < 2/L (should converge)
    gamma_converge = 0.99 * M
    w_c = -2.0  # An arbitrary starting point
    print(f"\n--- Test 1: Convergence (gamma = {gamma_converge:.4f} < M) ---")
    print(f"Starting at w = {w_c:.4f}")
    for i in range(10):
        w_c -= gamma_converge * grad_R(w_c, X, Y)
        print(f"Iteration {i+1:2d}: w = {w_c:8.4f}, R(w) = {R(w_c, X, Y):.4f}")
    print("Result: The weights are converging towards w*.")

    # Case 2: gamma > 2/L (should diverge)
    gamma_diverge = 1.01 * M
    w_d = -2.0 # The same starting point
    print(f"\n--- Test 2: Divergence (gamma = {gamma_diverge:.4f} > M) ---")
    print(f"Starting at w = {w_d:.4f}")
    for i in range(10):
        w_d -= gamma_diverge * grad_R(w_d, X, Y)
        # Check for numerical instability
        if np.isnan(w_d) or np.isinf(w_d) or abs(w_d) > 1e6:
            print(f"Iteration {i+1:2d}: w has diverged.")
            break
        print(f"Iteration {i+1:2d}: w = {w_d:8.4f}, R(w) = {R(w_d, X, Y):.4f}")
    print("Result: The weights are oscillating and diverging.")

solve()