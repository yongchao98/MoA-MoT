import numpy as np
from scipy.optimize import fsolve

def solve_micro_problem():
    """
    Solves for the elasticity of optimal search intensity q with respect to unemployment probability p.
    """
    # Initial parameter
    p = 0.2
    w = 20

    # This is the implicit function f(q, p) = 0 derived from the model's first-order conditions.
    # f(q, p) = exp(-2q) * (1.5 - 0.5p - q - p*q + p*q^2) - (1-p)*(1-q)
    def implicit_function(q, p_val):
        return np.exp(-2*q) * (1.5 - 0.5*p_val - q - p_val*q + p_val*q**2) - (1-p_val)*(1-q)

    # Step 1: Solve for the optimal q (q_star) for the given p = 0.2
    # We provide an initial guess of 0.5 for the root-finding algorithm.
    q_star = fsolve(implicit_function, 0.5, args=(p))[0]

    # Step 2: Calculate the partial derivatives of f(q, p) to find dq/dp.
    # ∂f/∂p = exp(-2q) * (q^2 - q - 0.5) + (1 - q)
    df_dp = np.exp(-2*q_star) * (q_star**2 - q_star - 0.5) + (1 - q_star)

    # ∂f/∂q = exp(-2q) * (-4 + 2q + 4pq - 2pq^2) + (1-p)
    df_dq = np.exp(-2*q_star) * (-4 + 2*q_star + 4*p*q_star - 2*p*q_star**2) + (1-p)

    # Step 3: Calculate dq/dp using the Implicit Function Theorem: dq/dp = - (∂f/∂p) / (∂f/∂q)
    dq_dp = -df_dp / df_dq

    # Step 4: Calculate the elasticity ε = (dq/dp) * (p / q)
    elasticity = dq_dp * (p / q_star)

    # Output the results step-by-step
    print(f"Given probability of unemployment, p = {p}")
    print(f"Calculated optimal search intensity, q* = {q_star:.5f}")
    print(f"The formula for elasticity is: ε = (dq/dp) * (p / q*)")
    print(f"The derivative dq/dp is calculated as: - (∂f/∂p) / (∂f/∂q) = -({df_dp:.5f}) / ({df_dq:.5f}) = {dq_dp:.5f}")
    print(f"Plugging in the numbers: ε = {dq_dp:.5f} * ({p} / {q_star:.5f})")
    print("\n---")
    print(f"Final Answer: The elasticity of optimal q with respect to p is {elasticity:.3f}")


solve_micro_problem()
<<<0.244>>>