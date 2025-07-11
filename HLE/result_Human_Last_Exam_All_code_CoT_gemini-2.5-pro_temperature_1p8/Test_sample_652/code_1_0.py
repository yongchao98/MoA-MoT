import numpy as np
from scipy.optimize import fsolve

def solve_and_calculate_elasticity():
    """
    This function solves for the optimal re-employment probability q and calculates
    its elasticity with respect to the unemployment probability p.
    """
    # Given parameter
    p = 0.2

    # Step 1: Define the equilibrium condition G(q, p) = 0.
    # This equation is the first-order condition of the government's problem.
    # G(q, p) = (1-q)e^(2q) - p(1-q)^2 - (1-p)(1.5-q) = 0
    G = lambda q: (1-q)*np.exp(2*q) - p*(1-q)**2 - (1-p)*(1.5-q)

    # Step 2: Solve for the optimal q* numerically for p = 0.2
    # We use an initial guess of 0.15 based on manual exploration.
    q_star = fsolve(G, 0.15)[0]

    print(f"The analysis is performed at the given unemployment probability p = {p:.1f}.")
    print("The equilibrium condition for the optimal re-employment probability q is: (1-q)e^(2q) - p(1-q)^2 - (1-p)(1.5-q) = 0")
    print(f"Solving this equation numerically gives the optimal q* = {q_star:.6f}\n")

    # Step 3: Define the partial derivatives of G with respect to q and p.
    # dG/dp = -q^2 + q + 0.5
    # dG/dq = e^(2q)*(1-2q) + 2p(1-q) + (1-p)
    dG_dp = -q_star**2 + q_star + 0.5
    dG_dq = np.exp(2*q_star)*(1 - 2*q_star) + 2*p*(1-q_star) + (1-p)

    print("To find the elasticity, we use the Implicit Function Theorem: dq/dp = - (dG/dp) / (dG/dq).")
    print(f"Partial derivative dG/dp at (p, q*) = {dG_dp:.6f}")
    print(f"Partial derivative dG/dq at (p, q*) = {dG_dq:.6f}\n")


    # Step 4: Calculate dq/dp using the Implicit Function Theorem.
    dq_dp = -dG_dp / dG_dq
    print(f"The derivative dq/dp is calculated as: -({dG_dp:.4f}) / ({dG_dq:.4f}) = {dq_dp:.6f}\n")

    # Step 5: Calculate the elasticity E = (dq/dp) * (p / q*)
    elasticity = dq_dp * (p / q_star)

    print("Finally, we calculate the elasticity of q with respect to p:")
    print(f"Elasticity E = (dq/dp) * (p / q*)")
    # Print the final equation with the numbers plugged in
    print(f"Elasticity E = ({dq_dp:.3f}) * ({p:.1f} / {q_star:.3f}) = {elasticity:.3f}")

if __name__ == "__main__":
    solve_and_calculate_elasticity()
