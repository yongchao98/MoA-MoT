import numpy as np
from scipy.optimize import minimize

def sigma(t):
    """Numerically stable sigmoid function."""
    t = np.clip(t, -500, 500)
    return 1 / (1 + np.exp(-t))

def R(w, x, y):
    """
    Objective function R(w).
    R(w) = 1/N * sum(log(1 + exp(-y_i * w * x_i)))
    """
    z = y * w * x
    return np.mean(np.log(1 + np.exp(-z)))

def R_grad(w, x, y):
    """Gradient of R(w)."""
    z = y * w * x
    return np.mean(y * x * sigma(-z))

def R_hess(w, x, y):
    """Hessian (second derivative) of R(w)."""
    z = y * w * x
    return np.mean(x**2 * sigma(z) * sigma(-z))

# Create a sample 1D non-separable dataset
x_data = np.array([-2.0, -1.0, 1.0, 2.5])
y_data = np.array([-1.0, 1.0, 1.0, -1.0])
N = len(x_data)

# Find the optimal w_* that minimizes R(w) using a numerical optimizer
result = minimize(
    fun=R,
    x0=0.0,  # Initial guess for w
    args=(x_data, y_data),
    jac=R_grad,
    method='L-BFGS-B'
)
w_star = result.x[0]

# Calculate the uniform (global) smoothness constant L
# L = sup_w R''(w) = R''(0) = (1 / (4 * N)) * sum(x_i^2)
L = np.mean(x_data**2) / 4.0

# Calculate the smoothness at the optimum, lambda
# lambda = R''(w_*)
lambda_val = R_hess(w_star, x_data, y_data)

# The theory for L-smooth convex functions guarantees convergence for gamma < 2/L.
# Thus, the largest upper bound M is 2/L.
M = 2 / L

print("Analysis of Gradient Descent Convergence for 1D Logistic Regression")
print("-" * 60)
print(f"Using sample dataset x = {x_data}, y = {y_data}")
print(f"The optimal weight found is w* = {w_star:.4f}")
print(f"The global smoothness constant is L = sup R''(w) = {L:.4f}")
print(f"The smoothness at the optimum is lambda = R''(w*) = {lambda_val:.4f}")
print("-" * 60)
print("For gradient descent to be guaranteed to converge for any initialization,")
print("the learning rate gamma must be less than M = 2 / L.")
print("\nFinal Equation:")
print(f"M = 2 / L = 2 / {L:.4f} = {M:.4f}")
print("\nThis corresponds to answer choice C.")
