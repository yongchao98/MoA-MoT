import numpy as np
from scipy.optimize import minimize

# Step 1: Define a sample non-separable 1D dataset (x_i, y_i)
# We choose 4 points. The points (x=-1, y=1) and (x=1, y=-1) make the set non-separable by a simple threshold at 0.
X = np.array([-2.0, -1.0, 1.0, 3.0])
Y = np.array([-1.0, 1.0, -1.0, 1.0])
N = len(X)

print(f"Using a sample dataset of N={N} points.")
print(f"X = {X}")
print(f"Y = {Y}\n")

# Step 2: Define the logistic risk function R(w), its gradient, and its second derivative
def sigma(t):
    """Numerically stable sigmoid function."""
    # Clip t to avoid overflow in exp
    t = np.clip(t, -500, 500)
    return 1 / (1 + np.exp(-t))

def R(w):
    """Logistic risk function R(w)."""
    # Use log(sigma(t)) = -log(1+exp(-t)) for stability
    # But for illustration, we can keep it simple as the values won't be extreme.
    # We add a small epsilon to avoid log(0) if sigma evaluates to 0, though unlikely here.
    return -np.mean(np.log(sigma(Y * w * X) + 1e-10))

def R_prime(w):
    """Gradient of R(w). Using the form sigma(-t) = 1 - sigma(t) is better."""
    return -np.mean(sigma(-Y * w * X) * Y * X)

def R_double_prime(w):
    """Second derivative (Hessian) of R(w)."""
    s_prime = sigma(-Y * w * X) * (1 - sigma(-Y * w * X))
    return np.mean(s_prime * X**2)

# Step 3: Calculate the smoothness constant L
# L is the supremum (maximum value) of R''(w).
# The function sigma'(t) = sigma(t)(1-sigma(t)) has a maximum value of 1/4 at t=0.
# Therefore, R''(w) is maximum at w=0.
L = R_double_prime(0)
# Theoretical formula: L = (1 / (4 * N)) * sum(X_i^2)
L_theory = (1 / (4 * N)) * np.sum(X**2)

print("--- Calculating Smoothness Constants ---")
print(f"The uniform smoothness constant L is the maximum value of the second derivative R''(w).")
print(f"This occurs at w=0. L = R''(0).")
print(f"Calculated L = {L:.4f}")
# print(f"Theoretical L = {L_theory:.4f}") # This should match

# Step 4: Find the optimal w* and calculate lambda
# We find w* by minimizing R(w) using a numerical optimizer.
result = minimize(R, x0=0.0, jac=R_prime)
w_star = result.x[0]

# Lambda is the smoothness (second derivative) at the optimal point w*
lambda_val = R_double_prime(w_star)

print(f"\nThe optimal weight w* is found by minimizing R(w).")
print(f"w* = {w_star:.4f}")
print(f"The smoothness at the optimum is lambda = R''(w*).")
print(f"Calculated lambda = {lambda_val:.4f}\n")
print("Note that lambda <= L, as expected.")


# Step 5: Calculate the largest upper bound M for the learning rate
# The theory of gradient descent for L-smooth convex functions states that convergence
# is guaranteed for a learning rate gamma < 2/L.
# Thus, the largest upper bound M is 2/L.
M = 2 / L

print("\n--- Determining the Convergence Bound M ---")
print("For guaranteed convergence from any initialization, the learning rate gamma must be less than 2/L.")
print("The largest upper bound M is therefore given by the expression: 2/L.")

# Outputting each number in the final equation
print("\nFor our example:")
print(f"The final equation for M is 2 / L")
print(f"M = 2 / {L:.4f} = {M:.4f}")

# For comparison, let's calculate the values for other options
bound_A = 1 / L
bound_B = 2 / lambda_val

print("\nComparing with other potential bounds:")
print(f"A. 1/L = {bound_A:.4f}")
print(f"B. 2/lambda = {bound_B:.4f}")
print(f"C. 2/L = {M:.4f}")
print("If we were to use a learning rate from (2/L, 2/lambda), e.g., gamma=1.8, we would risk divergence if starting near w=0.")
