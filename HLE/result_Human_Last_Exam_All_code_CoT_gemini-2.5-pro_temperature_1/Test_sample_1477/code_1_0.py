import numpy as np
from scipy.optimize import minimize

def sigma(t):
    """Numerically stable sigmoid function."""
    return np.where(t >= 0, 1 / (1 + np.exp(-t)), np.exp(t) / (1 + np.exp(t)))

def R(w, x, y):
    """Calculates the logistic loss R(w)."""
    N = len(x)
    # The loss is -1/N * sum(log(sigma(y_i * w * x_i)))
    # A stable way to compute log(sigma(t)) is -log(1 + exp(-t))
    z = y * w * x
    log_likelihood = -np.log(1 + np.exp(-z))
    return -np.sum(log_likelihood) / N

def R_prime(w, x, y):
    """Calculates the gradient (first derivative) of R(w)."""
    N = len(x)
    z = y * w * x
    # R'(w) = -1/N * sum(y_i * x_i * sigma(-z))
    return -np.sum(y * x * sigma(-z)) / N

def R_double_prime(w, x, y):
    """Calculates the Hessian (second derivative) of R(w)."""
    N = len(x)
    z = y * w * x
    # R''(w) = 1/N * sum(x_i^2 * sigma(z) * sigma(-z))
    return np.sum(x**2 * sigma(z) * sigma(-z)) / N

def run_gradient_descent(w_init, x, y, gamma, iterations):
    """Runs gradient descent for a given number of iterations."""
    w = w_init
    print(f"\nRunning GD with learning rate gamma = {gamma:.4f}")
    print("Iteration | Loss")
    print("----------|-----------")
    for i in range(iterations):
        loss = R(w, x, y)
        if i % (iterations // 5) == 0:
            print(f"{i:^10}| {loss:.7f}")
        grad = R_prime(w, x, y)
        w = w - gamma * grad
    final_loss = R(w, x, y)
    print(f"{iterations:^10}| {final_loss:.7f}")
    print("-" * 27)


# 1. Define a non-separable dataset
X = np.array([-2.0, -1.0, 1.0, 3.0])
y = np.array([1.0, -1.0, 1.0, -1.0])
N = len(X)
w_init = 0.0

# 2. Calculate L (global smoothness constant)
# L = sup(R''(w)) = sum(x_i^2) / (4*N)
L = np.sum(X**2) / (4 * N)

# 3. Find the optimal w* to calculate lambda
# Use a numerical optimizer to find the minimum of R(w)
res = minimize(fun=R, x0=w_init, args=(X, y), jac=R_prime, method='L-BFGS-B')
w_star = res.x[0]

# 4. Calculate lambda (smoothness at the optimum)
lambda_val = R_double_prime(w_star, X, y)

# 5. Define the convergence boundary M
M = 2 / L

# 6. Print the theoretical values
print("--- Theoretical Analysis ---")
print(f"Dataset X: {X}")
print(f"Dataset y: {y}")
print(f"L (Global smoothness)      = sum(x_i^2)/(4N) = {np.sum(X**2)}/(4*{N}) = {L:.4f}")
print(f"Optimal w* found           = {w_star:.4f}")
print(f"lambda (Smoothness at w*)  = R''(w*) = {lambda_val:.4f}")
print("\n--- Final Equation for Convergence Bound ---")
print(f"M = 2 / L")
print(f"M = 2 / {L:.4f} = {M:.4f}")
print("-" * 36)


# 7. Demonstrate convergence vs. divergence
print("\n--- Gradient Descent Demonstration ---")

# Experiment 1: gamma < M (should converge)
gamma_safe = 0.99 * M
run_gradient_descent(w_init, X, y, gamma_safe, iterations=20)

# Experiment 2: gamma > M (should diverge or oscillate)
gamma_unsafe = 1.01 * M
run_gradient_descent(w_init, X, y, gamma_unsafe, iterations=20)
