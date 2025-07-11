import numpy as np
from scipy.optimize import root_scalar

# Define the value of p
p = 0.2

# Define the implicit function H(q, p) = 0 derived from the optimization problem.
# (1-q)exp(2q) - 1.5 + q - (p/(1-p))(1-q)^2 = 0
def implicit_H(q, p_val):
    """Implicit function H(q, p) that defines the optimal q."""
    p_term = p_val / (1 - p_val)
    return (1 - q) * np.exp(2 * q) - 1.5 + q - p_term * (1 - q)**2

# --- Step 1: Solve for optimal q* at p = 0.2 ---
# We use a numerical root finder to solve H(q, 0.2) = 0.
# The solution must be in (0,1). We can bracket the root between 0.3 and 0.4 from initial exploration.
sol = root_scalar(lambda q: implicit_H(q, p), bracket=[0.3, 0.4], method='brentq')
q_star = sol.root

# --- Step 2: Calculate the derivative dq/dp ---
# This is done using the Implicit Function Theorem: dq/dp = - (dH/dp) / (dH/dq)

# Define the partial derivative of H with respect to q
def dH_dq(q, p_val):
    """Partial derivative of H with respect to q."""
    p_term = p_val / (1 - p_val)
    return (1 - 2 * q) * np.exp(2 * q) + 1 + 2 * p_term * (1 - q)

# Define the partial derivative of H with respect to p
def dH_dp(q, p_val):
    """Partial derivative of H with respect to p."""
    return -((1 - q)**2) / ((1 - p_val)**2)

# Evaluate the partial derivatives at (q*, p)
dH_dp_val = dH_dp(q_star, p)
dH_dq_val = dH_dq(q_star, p)

# Calculate dq/dp
dq_dp = -dH_dp_val / dH_dq_val

# --- Step 3: Calculate the elasticity ---
# Elasticity = (dq/dp) * (p / q*)
elasticity = dq_dp * (p / q_star)

# --- Step 4: Print the results ---
print(f"For a probability of unemployment p = {p:.1f}:")
print(f"The optimal re-employment probability is q* = {q_star:.4f}")
print(f"The derivative of q with respect to p is dq/dp = {dq_dp:.4f}")
print("\nThe elasticity is calculated as (dq/dp) * (p / q*):")
print(f"Final Equation: {dq_dp:.4f} * ({p:.1f} / {q_star:.4f}) = {elasticity:.3f}")
