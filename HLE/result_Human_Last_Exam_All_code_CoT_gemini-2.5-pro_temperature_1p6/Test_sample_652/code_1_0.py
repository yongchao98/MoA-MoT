import numpy as np
from scipy.optimize import fsolve

def solve_elasticity():
    """
    This function solves the microeconomic model to find the elasticity of optimal
    job search intensity (q) with respect to the probability of unemployment (p).
    """
    # Step 1: Define model parameters
    w = 20.0
    p = 0.2

    # Step 2: Define the implicit function H(q, p) = 0, which characterizes the optimal q.
    # This equation is derived from the first-order conditions of the government's
    # utility maximization problem, subject to the worker's incentive compatibility
    # and the government's budget constraint.
    def H(q, p_val=p):
        """Implicit function H(q, p) = 0."""
        term1 = np.exp(2 * q)
        term2 = - p_val * (1 - q) / (1 - p_val)
        term3 = - (1.5 - q) / (1 - q)
        return term1 + term2 + term3

    # Step 3: Solve for the optimal q at p = 0.2 using a numerical solver.
    # q must be in (0, 1), so an initial guess of 0.5 is reasonable.
    q_initial_guess = 0.5
    optimal_q = fsolve(H, q_initial_guess, args=(p,))[0]

    # Step 4: Define the partial derivatives of H with respect to p and q.
    def dH_dp(q, p_val=p):
        """Partial derivative of H with respect to p."""
        # The derivative of -p*(1-q)/(1-p) w.r.t p is -(1-q)/(1-p)^2
        return -(1 - q) / ((1 - p_val)**2)

    def dH_dq(q, p_val=p):
        """Partial derivative of H with respect to q."""
        # Derivative of exp(2q) is 2*exp(2q)
        # Derivative of -p*(1-q)/(1-p) is p/(1-p)
        # Derivative of -(1.5-q)/(1-q) is -0.5/(1-q)^2
        term1 = 2 * np.exp(2 * q)
        term2 = p_val / (1 - p_val)
        term3 = -0.5 / ((1 - q)**2)
        return term1 + term2 + term3

    # Step 5: Calculate dq/dp using the Implicit Function Theorem.
    dH_dp_val = dH_dp(optimal_q, p)
    dH_dq_val = dH_dq(optimal_q, p)
    dq_dp = -dH_dp_val / dH_dq_val

    # Step 6: Calculate the elasticity of q with respect to p.
    elasticity = (dq_dp * p) / optimal_q
    
    print("Step-by-step calculation:")
    print(f"Given probability of unemployment, p = {p}")
    print(f"Optimal job search intensity, q* = {optimal_q:.4f}")
    print(f"Partial derivative ∂H/∂p at (q*, p) = {dH_dp_val:.4f}")
    print(f"Partial derivative ∂H/∂q at (q*, p) = {dH_dq_val:.4f}")
    print(f"Derivative dq/dp = -({dH_dp_val:.4f}) / ({dH_dq_val:.4f}) = {dq_dp:.4f}")
    print(f"Elasticity E = ({dq_dp:.4f} * {p}) / {optimal_q:.4f}")
    print(f"The elasticity of optimal q with respect to p is: {elasticity:.3f}")
    
    return elasticity

# Run the calculation and print the final answer in the required format
result = solve_elasticity()
# The final result is presented below
print(f"<<<{result:.3f}>>>")