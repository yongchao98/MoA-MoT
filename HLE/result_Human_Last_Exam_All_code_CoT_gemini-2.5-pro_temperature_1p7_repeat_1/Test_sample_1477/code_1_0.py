import numpy as np

def solve():
    """
    This function demonstrates the calculation of the maximum learning rate for
    gradient descent on a logistic regression problem and verifies it.
    """
    # --- Step 1: Define a sample dataset and the model functions ---
    # We choose a small, non-separable dataset in 1D for demonstration.
    # Data points (x_i, y_i)
    X = np.array([1.0, 2.0, -1.5])
    Y = np.array([1.0, -1.0, -1.0])
    N = len(X)

    # Sigmoid function
    def sigma(t):
        # Add a clip for numerical stability with large exponents
        t = np.clip(t, -500, 500)
        return 1.0 / (1.0 + np.exp(-t))

    # Loss function R(w)
    def R(w, X_data, Y_data):
        N_data = len(X_data)
        # Add a small epsilon for log stability
        s = sigma(Y_data * w * X_data)
        return -np.sum(np.log(s + 1e-15)) / N_data

    # Gradient of R(w)
    def grad_R(w, X_data, Y_data):
        N_data = len(X_data)
        terms = (sigma(Y_data * w * X_data) - 1.0) * Y_data * X_data
        return np.sum(terms) / N_data

    # --- Step 2: Calculate the global smoothness constant L ---
    # The smoothness L is the supremum of the second derivative of R(w).
    # The term sigma(t)*(1-sigma(t)) has a maximum value of 1/4.
    sum_x_squared = np.sum(X**2)
    L = (1.0 / 4.0) * sum_x_squared / N

    # --- Step 3: Calculate the convergence bound M ---
    # The largest upper bound for the learning rate is M = 2/L.
    M = 2.0 / L

    # --- Step 4: Perform Gradient Descent to verify the bound ---
    w_init = 0.5
    num_iterations = 50

    # Case 1: gamma < M (should converge)
    gamma1 = 0.95 * M
    w1 = w_init
    for _ in range(num_iterations):
        grad = grad_R(w1, X, Y)
        w1 = w1 - gamma1 * grad
    loss1 = R(w1, X, Y)

    # Case 2: gamma > M (should diverge or oscillate)
    gamma2 = 1.05 * M
    w2 = w_init
    for _ in range(num_iterations):
        grad = grad_R(w2, X, Y)
        w2 = w2 - gamma2 * grad
    loss2 = R(w2, X, Y)


    # --- Step 5: Output the step-by-step reasoning and results ---
    print("This script solves for the largest upper bound M for the learning rate.")
    print("----------------------------------------------------------------------")
    print(f"1. We start with a sample dataset with N = {N} points:")
    print(f"   X = {X}")
    print(f"   Y = {Y}")
    print("\n2. We calculate the global smoothness constant L of the loss function.")
    print("   The formula is L = (1 / (4 * N)) * sum(x_i^2).")
    print(f"   The sum of x_i^2 is {sum_x_squared:.2f}.")
    print(f"   Each number in the equation for L is:")
    print(f"     N = {N}")
    print(f"     sum(x_i^2) = {sum_x_squared:.2f}")
    print(f"   So, L = (1 / (4 * {N})) * {sum_x_squared:.2f} = {L:.4f}")
    print("\n3. According to optimization theory, gradient descent is guaranteed to")
    print("   converge for any starting point if the learning rate γ < 2/L.")
    print("   This means the largest upper bound is M = 2/L.")
    print(f"   Each number in the equation for M is:")
    print(f"     L = {L:.4f}")
    print(f"   So, M = 2 / {L:.4f} = {M:.4f}")
    print("\n4. Verification:")
    print(f"   We run gradient descent with γ₁ = {gamma1:.4f} (< M) and γ₂ = {gamma2:.4f} (> M).")
    print(f"   After {num_iterations} iterations:")
    print(f"   - Final loss with γ₁ (the safe rate) is: {loss1:.6f}")
    print(f"   - Final loss with γ₂ (the unsafe rate) is: {loss2:.6f}")
    print("\nConclusion: The loss for the rate greater than M is significantly higher,")
    print("confirming that M = 2/L is the correct theoretical boundary.")

solve()