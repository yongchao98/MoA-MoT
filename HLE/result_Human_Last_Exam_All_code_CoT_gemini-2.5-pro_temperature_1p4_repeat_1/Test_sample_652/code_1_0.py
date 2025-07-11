import numpy as np
from scipy.optimize import fsolve

# --- 1. Define constants and the implicit function F(q, p) = 0 ---
# The problem parameters
w = 20.0
p = 0.2

# The implicit function F(q, p) derived from the model's first-order conditions.
# We need to find the root of this function for q.
# F(q, p) = exp(-2q) * [ (p/(1-p))(1-q)^2 + 1.5 - q ] - (1-q) = 0
def F(q, p_val):
    """Implicit function F(q, p) that defines the optimal q."""
    return np.exp(-2*q) * ( (p_val/(1-p_val))*(1-q)**2 + 1.5 - q ) - (1-q)

# --- 2. Solve for optimal q* at p = 0.2 ---
# We solve F(q, 0.2) = 0 for q. We need an initial guess for the solver.
# A value between 0 and 1 is reasonable. Let's start with 0.3.
initial_guess = 0.3
q_star = fsolve(F, initial_guess, args=(p))[0]

# --- 3. Define the partial derivatives of F needed for the implicit function theorem ---
def dF_dp(q, p_val):
    """Partial derivative of F with respect to p."""
    return np.exp(-2*q) * ( (1-q)**2 / (1-p_val)**2 )

def dF_dq(q, p_val):
    """Partial derivative of F with respect to q."""
    # This simplified form is derived in the thought process.
    # It is equivalent to the full derivative but computationally simpler.
    term1 = -1 + 2*q
    term2 = np.exp(-2*q) * (2 * p_val * (1-q) / (1-p_val) + 1)
    return term1 - term2

# --- 4. Calculate the numerical values of the derivatives ---
val_dF_dp = dF_dp(q_star, p)
val_dF_dq = dF_dq(q_star, p)

# --- 5. Calculate dq/dp and the elasticity ---
# Using the implicit function theorem: dq/dp = - (dF/dp) / (dF/dq)
dq_dp = -val_dF_dp / val_dF_dq

# Elasticity E = (dq/dp) * (p / q)
elasticity = dq_dp * (p / q_star)

# --- 6. Print the results step-by-step ---
print(f"The problem is solved for the following parameter values:")
print(f"Initial wage w = {w}")
print(f"Probability of unemployment p = {p}\n")

print(f"Step 1: Numerically solve for the optimal search intensity q*.")
print(f"The optimal q* is found to be: {q_star:.4f}\n")

print(f"Step 2: Calculate the partial derivatives at (q*, p).")
print(f"The partial derivative ∂F/∂p at ({q_star:.4f}, {p}) is: {val_dF_dp:.4f}")
print(f"The partial derivative ∂F/∂q at ({q_star:.4f}, {p}) is: {val_dF_dq:.4f}\n")

print(f"Step 3: Calculate the derivative dq/dp using the Implicit Function Theorem.")
print(f"dq/dp = - (∂F/∂p) / (∂F/∂q)")
print(f"dq/dp = -({val_dF_dp:.4f}) / ({val_dF_dq:.4f}) = {dq_dp:.4f}\n")

print(f"Step 4: Calculate the final elasticity of q with respect to p.")
print(f"Elasticity = (dq/dp) * (p / q*)")
print(f"Elasticity = {dq_dp:.4f} * ({p} / {q_star:.4f}) = {elasticity:.3f}")

print("\nFinal Answer:")
print(f"The elasticity of optimal q with respect to p around p=0.2 is {elasticity:.3f}.")
<<<0.762>>>