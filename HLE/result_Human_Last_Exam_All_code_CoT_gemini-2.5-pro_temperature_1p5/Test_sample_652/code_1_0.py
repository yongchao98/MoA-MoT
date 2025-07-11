import numpy as np
from scipy.optimize import fsolve

# --- Step 1: Define the model parameters and the implicit function ---

# Given probability of unemployment
p_val = 0.2

# The implicit function F(q, p) = 0 derived from the first-order condition of the government's problem.
# (1-q) - exp(-2q) * [ (1.5 - q) + p*(1-q)^2 ] = 0
def F(q, p):
    """Implicit function for the optimal q."""
    term1 = 1 - q
    term2 = np.exp(-2 * q)
    term3 = (1.5 - q) + p * (1 - q)**2
    return term1 - term2 * term3

# --- Step 2: Solve for the optimal q* at p = 0.2 ---

# We need an initial guess for q. Since q is a probability, it's between 0 and 1.
# Let's try 0.5 as a starting point.
q_initial_guess = 0.5
# Use fsolve to find the root of F(q, p_val) = 0
q_star = fsolve(F, q_initial_guess, args=(p_val,))[0]


# --- Step 3: Define the partial derivatives for the Implicit Function Theorem ---

def dF_dp(q, p):
    """Partial derivative of F with respect to p."""
    return -np.exp(-2 * q) * (1 - q)**2

def dF_dq(q, p):
    """Partial derivative of F with respect to q."""
    exp_term = np.exp(-2 * q)
    # The derivative can be calculated as:
    # dF/dq = -1 + 2*exp(-2q)*[(1.5 - q) + p*(1-q)**2] + exp(-2q)*[1 + 2*p*(1-q)]
    term_A = (1.5 - q) + p * (1 - q)**2
    term_B = 1 + 2 * p * (1 - q)
    return -1 + 2 * exp_term * term_A + exp_term * term_B

# --- Step 4: Calculate dq/dp and the elasticity ---

# Evaluate the partial derivatives at the optimal point (q_star, p_val)
dFdp_val = dF_dp(q_star, p_val)
dFdq_val = dF_dq(q_star, p_val)

# Calculate dq/dp using the Implicit Function Theorem: dq/dp = - (dF/dp) / (dF/dq)
dq_dp = -dFdp_val / dFdq_val

# Calculate the elasticity E = (dq/dp) * (p / q*)
elasticity = dq_dp * (p_val / q_star)

# --- Step 5: Print the results ---

print("--- Calculation Steps ---")
print(f"1. Optimal search intensity (q*) for p={p_val}: {q_star:.5f}")
print(f"2. Sensitivity of q* to p (dq/dp): {dq_dp:.5f}")
print("\n--- Final Elasticity Calculation ---")
print(f"Elasticity = (dq/dp) * (p / q*)")
print(f"           = {dq_dp:.5f} * ({p_val} / {q_star:.5f})")
print(f"\nThe elasticity of optimal q with respect to p is: {elasticity:.3f}")

print(f"\n<<<{elasticity:.3f}>>>")