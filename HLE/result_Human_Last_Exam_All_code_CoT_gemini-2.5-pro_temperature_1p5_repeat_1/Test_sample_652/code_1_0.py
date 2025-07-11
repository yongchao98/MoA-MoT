import numpy as np
from scipy.optimize import fsolve

def solve_and_calculate_elasticity():
    """
    This function solves the microeconomic model to find the elasticity of optimal 
    job search intensity (q) with respect to the probability of unemployment (p).
    """
    # The value of p around which we calculate the elasticity
    p_val = 0.2

    # The equilibrium condition derived from the model is H(q, p) = 0.
    # H(q, p) = exp(2*q) - (1.5 - q)/(1 - q) - (p/(1-p))*(1-q)
    def H(q, p):
        # The search intensity q must be a probability between 0 and 1.
        if q <= 0 or q >= 1:
            return np.inf
        
        term1 = np.exp(2 * q)
        term2 = (1.5 - q) / (1 - q)
        term3 = (p / (1 - p)) * (1 - q)
        return term1 - term2 - term3

    # We need to find the optimal q for p=0.2 by finding the root of H(q, 0.2) = 0.
    # We use a numerical solver, fsolve, with an initial guess for q.
    q_initial_guess = 0.3
    q_star = fsolve(H, q_initial_guess, args=(p_val,))[0]

    # To calculate the elasticity, we use the Implicit Function Theorem:
    # E = (dq/dp) * (p/q), where dq/dp = - (∂H/∂p) / (∂H/∂q).
    # We define the functions for the partial derivatives of H.

    # Partial derivative of H with respect to p
    def dH_dp(q, p):
        return -(1 - q) / ((1 - p)**2)

    # Partial derivative of H with respect to q
    def dH_dq(q, p):
        term1 = 2 * np.exp(2 * q)
        term2 = 0.5 / ((1 - q)**2)
        term3 = p / (1 - p)
        return term1 - term2 + term3

    # Now, we evaluate these derivatives at the equilibrium point (q_star, p_val).
    dH_dp_val = dH_dp(q_star, p_val)
    dH_dq_val = dH_dq(q_star, p_val)
    
    # Calculate dq/dp
    dq_dp = -dH_dp_val / dH_dq_val
    
    # Finally, calculate the elasticity.
    elasticity = dq_dp * (p_val / q_star)

    # --- Output the results step-by-step ---
    print(f"The analysis leads to an equilibrium condition H(q, p) = 0.")
    print(f"For p = {p_val}, we solve this equation numerically.")
    print(f"The optimal search intensity is q* = {q_star:.5f}\n")
    
    print("To find the elasticity E = (dq/dp) * (p/q), we calculate the components:")
    print("1. The partial derivatives of H(q, p) at the equilibrium point:")
    print(f"   ∂H/∂p = {dH_dp_val:.5f}")
    print(f"   ∂H/∂q = {dH_dq_val:.5f}\n")
    
    print("2. The derivative dq/dp using the Implicit Function Theorem:")
    print(f"   dq/dp = -({dH_dp_val:.5f}) / ({dH_dq_val:.5f}) = {dq_dp:.5f}\n")
    
    print("3. The final elasticity calculation:")
    print(f"   Elasticity = dq/dp * (p/q*)")
    print(f"   Elasticity = {dq_dp:.5f} * ({p_val} / {q_star:.5f})")
    print(f"   Elasticity = {elasticity:.3f}\n")

# Run the calculation and print the results
solve_and_calculate_elasticity()