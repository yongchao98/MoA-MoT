import numpy as np
from scipy.optimize import fsolve

# --- 1. Define problem parameters ---
p = 0.2

# --- 2. Find the optimal q for p=0.2 ---
# The equilibrium condition that implicitly defines q as a function of p is F(q,p) = 0.
# F(q, p) = exp(-2q)*[1.5 - q + (p/(1-p))*(1-q)^2] - (1-q)
# For p=0.2, this simplifies to: exp(-2q)*(q^2 - 6q + 7) - 4*(1-q) = 0
def equilibrium_equation(q, p_val):
    """
    The first-order condition F(q, p) = 0 that defines the optimal q.
    """
    K = p_val / (1 - p_val)
    # Using the original form to be general
    term1 = np.exp(-2*q) * (1.5 - q + K * (1-q)**2)
    term2 = -(1-q)
    return term1 + term2

# We need to find the root of this equation for p=0.2.
# Initial guess for q (must be between 0 and 1)
q_initial_guess = 0.3
# Use a numerical solver to find the value of q
q_optimal = fsolve(lambda q: equilibrium_equation(q, p), q_initial_guess)[0]

# --- 3. Calculate the derivative dq/dp using the Implicit Function Theorem ---
# dq/dp = - (dF/dp) / (dF/dq)

def dF_dp(q, p_val):
    """Calculates the partial derivative of F with respect to p."""
    # Derivative of p/(1-p) is 1/(1-p)^2
    dK_dp = 1 / (1 - p_val)**2
    return dK_dp * np.exp(-2*q) * (1-q)**2

def dF_dq(q, p_val):
    """Calculates the partial derivative of F with respect to q."""
    K = p_val / (1 - p_val)
    # Derivative of exp(-2q)*(1.5-q) + 1
    part1 = np.exp(-2*q) * (2*q - 4) + 1
    # Derivative of K*exp(-2q)*(1-q)^2
    part2 = K * (-2 * np.exp(-2*q) * (1-q)**2 + np.exp(-2*q) * 2 * (1-q) * (-1))
    return part1 + part2

# Evaluate the partial derivatives at the optimal q and given p
dFdp_val = dF_dp(q_optimal, p)
dFdq_val = dF_dq(q_optimal, p)

# Calculate dq/dp
dq_dp = -dFdp_val / dFdq_val

# --- 4. Calculate the elasticity ---
# Elasticity E = (dq/dp) * (p/q)
elasticity = dq_dp * (p / q_optimal)

# --- 5. Print the results ---
print("This script calculates the elasticity of the optimal job search intensity (q) with respect to the probability of unemployment (p).")
print("\nThe final elasticity is calculated using the formula: E = (dq/dp) * (p / q)")
print("\nFirst, we solve for the optimal 'q' at p = 0.2:")
print(f"q = {q_optimal:.5f}")
print("\nNext, we calculate the derivative dq/dp using the implicit function theorem:")
print(f"dq/dp = {dq_dp:.5f}")
print("\nFinally, we plug these numbers into the elasticity formula:")
print(f"Elasticity = {dq_dp:.5f} * ({p:.1f} / {q_optimal:.5f})")
print(f"\nThe calculated elasticity is: {elasticity:.3f}")
