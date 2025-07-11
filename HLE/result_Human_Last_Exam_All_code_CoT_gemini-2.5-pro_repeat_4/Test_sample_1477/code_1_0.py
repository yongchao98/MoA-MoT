import numpy as np

# This script demonstrates the convergence condition for gradient descent
# on a logistic regression problem.

# Step 1: Define a non-separable dataset and the corresponding loss function R(w)
# and its derivatives.
# Dataset: (x1, y1) = (1, 1), (x2, y2) = (2, -1)
# N = 2
x = np.array([1., 2.])
y = np.array([1., -1.])

# Loss function: R(w) = 0.5 * [log(1 + exp(-w)) + log(1 + exp(2w))]
# Note: We work with the sum of losses, so N=1 in the code implementation.
# This scales R and its derivatives by a factor of N=2, which is fine as
# the learning rate condition gamma < 2/L will scale accordingly and the logic holds.
# For consistency with the derivation, we can divide by N=2.
N = 2.0

def grad_R(w):
    """Computes the gradient of R(w)."""
    term1_grad = -y[0] * x[0] / (1 + np.exp(y[0] * x[0] * w))
    term2_grad = -y[1] * x[1] / (1 + np.exp(y[1] * x[1] * w))
    return (term1_grad + term2_grad) / N

def hessian_R(w):
    """Computes the second derivative (Hessian in 1D) of R(w)."""
    # sigma(t)*(1-sigma(t)) = exp(-t) / (1+exp(-t))^2 = exp(t) / (1+exp(t))^2
    term1_hess = (x[0]**2) * np.exp(y[0] * x[0] * w) / (1 + np.exp(y[0] * x[0] * w))**2
    term2_hess = (x[1]**2) * np.exp(y[1] * x[1] * w) / (1 + np.exp(y[1] * x[1] * w))**2
    return (term1_hess + term2_hess) / N

# Step 2: Calculate the global smoothness L and the critical learning rate M.
# L is the maximum value of the second derivative, which occurs at w=0.
L = hessian_R(0)
M = 2 / L

print("--- Theoretical Analysis ---")
print(f"The smoothness constant L = sup|R''(w)| is: {L:.4f}")
print("The largest upper bound M for the learning rate is 2/L.")
print(f"M = 2 / {L:.4f} = {M:.4f}")
print("Any learning rate gamma < M should converge.")
print("Any learning rate gamma > M may diverge.")
print("-" * 28 + "\n")


# Step 3: Run Gradient Descent with different learning rates.
def gradient_descent(w_init, gamma, num_iterations):
    """Performs gradient descent."""
    w = w_init
    w_history = [w]
    print(f"Running GD with gamma = {gamma:.4f}, starting from w_init = {w_init:.2f}")
    for i in range(num_iterations):
        gradient = grad_R(w)
        w = w - gamma * gradient
        w_history.append(w)
        if i < 10 or (i > 90 and i < 100): # Print some early and late steps
            print(f"  Iter {i+1:3d}: w = {w:.6f}")
        if abs(w) > 1e6:
            print("  Divergence detected!")
            break
    return w_history

# Case 1: gamma < M (should converge)
gamma_converge = M - 0.1
w_init = 0.1
print("--- Case 1: gamma < M (Convergence) ---")
gradient_descent(w_init, gamma_converge, 100)
print("-" * 40 + "\n")

# Case 2: gamma > M (should diverge)
gamma_diverge = M + 0.1
print("--- Case 2: gamma > M (Divergence) ---")
gradient_descent(w_init, gamma_diverge, 100)
print("-" * 40 + "\n")
