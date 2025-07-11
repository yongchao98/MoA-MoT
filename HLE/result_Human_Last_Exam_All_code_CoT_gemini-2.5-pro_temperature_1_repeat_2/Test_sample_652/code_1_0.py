import numpy as np
from scipy.optimize import fsolve

# --- Step 1 & 2: Define parameters and the equilibrium equation ---
# The problem parameters
p_val = 0.2
w = 20

# The equilibrium is defined by the implicit function G(q, p) = 0, derived from the model's FOCs.
# G(q, p) = exp(-2q) * (1 - p*q) - (1-q)
def G(q, p):
    """Implicit function G(q, p) = 0 that defines the equilibrium q."""
    return np.exp(-2 * q) * (1 - p * q) - (1 - q)

# --- Step 3: Solve for optimal q at p=0.2 ---
# We need an initial guess for the root finder. A value between 0 and 1 is appropriate.
initial_guess_q = 0.8
# Use fsolve to find the root q_val that satisfies G(q_val, p_val) = 0
q_val = fsolve(G, initial_guess_q, args=(p_val,))[0]


# --- Step 4: Calculate dq/dp using the Implicit Function Theorem ---
# dq/dp = - (dG/dp) / (dG/dq)

# Partial derivative of G with respect to p
# dG/dp = d/dp [exp(-2q)*(1-pq) - (1-q)] = exp(-2q)*(-q)
dG_dp = -q_val * np.exp(-2 * q_val)

# Partial derivative of G with respect to q
# dG/dq = d/dq [exp(-2q)*(1-pq) - (1-q)] = -2*exp(-2q)*(1-pq) + exp(-2q)*(-p) - (-1)
# dG/dq = exp(-2q) * (-2(1-pq) - p) + 1
dG_dq = np.exp(-2 * q_val) * (-2 * (1 - p_val * q_val) - p_val) + 1

# Calculate dq/dp
dq_dp = -dG_dp / dG_dq

# --- Step 5: Calculate the elasticity ---
# Elasticity E = (dq/dp) * (p/q)
elasticity = dq_dp * (p_val / q_val)

# --- Step 6: Print the results ---
# The final equation is E = (dq/dp) * (p/q)
print("To calculate the elasticity, we first find the values for each component:")
print(f"Given p = {p_val:.3f}")
print(f"The optimal search intensity is q = {q_val:.3f}")
print(f"The change of q with respect to p is dq/dp = {dq_dp:.3f}")
print("\nPlugging these into the elasticity formula:")
print(f"Elasticity = (dq/dp) * (p / q)")
print(f"           = {dq_dp:.3f} * ({p_val:.3f} / {q_val:.3f})")
print(f"           = {elasticity:.3f}")

print(f"<<<{elasticity:.3f}>>>")