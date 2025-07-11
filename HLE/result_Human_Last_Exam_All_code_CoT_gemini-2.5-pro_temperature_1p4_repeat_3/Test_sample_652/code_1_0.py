import numpy as np
from scipy.optimize import fsolve

def solve_and_calculate_elasticity():
    """
    This function solves for the optimal job search intensity q and calculates its
    elasticity with respect to the unemployment probability p.
    """
    # Step 1: Define the parameters from the problem.
    # We are calculating the elasticity around p = 0.2.
    p = 0.2

    # Step 2: Define the implicit function G(q, p) = 0.
    # The optimal q as a function of p, q*(p), is the root of this equation.
    # This equation is derived from the first-order condition of the government's
    # utility maximization problem.
    # The equation is: e^{-2q}[ (1.5-q) + p(1-q)^2/(1-p) ] - (1-q) = 0
    def G(q, p_val):
        # We handle cases where q is out of the valid (0,1) range.
        if q >= 1 or q <= 0:
            return np.inf
        
        term1 = 1.5 - q
        term2 = p_val * (1 - q)**2 / (1 - p_val)
        return np.exp(-2 * q) * (term1 + term2) - (1 - q)

    # Step 3: Solve for the optimal q (q_star) at p=0.2.
    # We use a numerical solver (fsolve) with an initial guess of 0.3.
    # The function passed to fsolve must be a function of only one variable, q.
    q_star = fsolve(lambda q: G(q, p), 0.3)[0]

    # Step 4: Calculate the partial derivatives for implicit differentiation.
    # The formula is dq/dp = - (dG/dp) / (dG/dq).

    # Partial derivative of G with respect to p:
    # dG/dp = exp(-2q) * (1-q)^2 / (1-p)^2
    dG_dp = np.exp(-2 * q_star) * (1 - q_star)**2 / (1 - p)**2

    # Partial derivative of G with respect to q:
    # dG/dq = 2*q - 1 - exp(-2*q) * (1 + 2*p*(1-q)/(1-p))
    dG_dq = 2 * q_star - 1 - np.exp(-2 * q_star) * (1 + 2 * p * (1 - q_star) / (1 - p))

    # Step 5: Calculate the derivative dq*/dp.
    dq_dp = -dG_dp / dG_dq

    # Step 6: Calculate the elasticity of q with respect to p.
    # Elasticity = (dq/dp) * (p / q)
    elasticity = dq_dp * (p / q_star)

    # Step 7: Print the components of the final equation and the result.
    print("The final elasticity E is calculated using the equation E = (dq*/dp) * (p / q*).")
    print("The values used are:")
    print(f"p = {p}")
    print(f"Optimal q* = {q_star:.5f}")
    print(f"Derivative dq*/dp = {dq_dp:.5f}")
    print(f"\nFinal calculated elasticity:")
    print(f"{elasticity:.3f}")

solve_and_calculate_elasticity()