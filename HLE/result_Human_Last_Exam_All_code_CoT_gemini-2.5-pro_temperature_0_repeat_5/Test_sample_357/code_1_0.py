import math

def solve_contract_problem():
    """
    This function calculates the separating equilibrium contracts for Chelsea's handyman problem.
    """
    # Step 1: Determine the optimal contract for the low type (L).
    # Chelsea's profit simplifies to (1/6) * (v_L - v_L^2) when IR constraints are binding.
    # To maximize profit, we maximize the function f(v_L) = v_L - v_L^2.
    # The maximum of this quadratic function occurs at v_L = -b / (2a), where a=-1 and b=1.
    v_L = -1 / (2 * -1)

    # The corresponding price p_L is determined by the binding IR_L constraint: p_L = v_L^2.
    p_L = v_L**2

    # Step 2: Determine the optimal contract for the high type (H).
    # Chelsea's profit is independent of v_H. We assume she chooses the highest possible v_H
    # to maximize the value added to her home, subject to constraints.
    # The constraints on v_H are v_H <= 30 (productivity) and v_H >= 1 (from IC_L).
    # Therefore, the optimal v_H is the maximum value in the interval [1, 30].
    v_H = 30.0

    # The corresponding price p_H is determined by the binding IR_H constraint: p_H = v_H.
    p_H = v_H

    # Step 3: Print the resulting pair of contracts.
    print("The separating equilibrium consists of the following pair of contracts:")
    print(f"Contract for Low Type (L): (v_L = {v_L}, p_L = {p_L})")
    print(f"Contract for High Type (H): (v_H = {v_H}, p_H = {p_H})")

solve_contract_problem()