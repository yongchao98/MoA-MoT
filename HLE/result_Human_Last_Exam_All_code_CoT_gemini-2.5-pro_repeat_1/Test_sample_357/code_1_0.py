import math

def solve_contract_problem():
    """
    This function solves for the separating equilibrium contracts (v_L, p_L) and (v_H, p_H).

    The process involves setting up and solving the principal's (Chelsea's) optimization problem
    subject to the participation and incentive compatibility constraints for both handyman types.
    """

    # Step 1: Define constraints and simplify
    # Based on economic theory, the low type's participation constraint (IR_L) and the high
    # type's incentive compatibility constraint (IC_H) will be binding.
    #
    # IR_L (binding): p_L - v_L^2 = 0  =>  p_L = v_L^2
    # IC_H (binding): p_H - v_H = p_L - v_L  =>  p_H = v_H + p_L - v_L

    # Step 2: Solve for the optimal v_L for Chelsea
    # Chelsea's expected profit simplifies to E[π] = v_L - v_L^2.
    # To maximize this, she must consider the constraints on v_L.
    # The productivity constraint for the low type is: v_L <= 20 - (1/3)*v_L  =>  v_L <= 15.
    # A second constraint arises from checking the non-binding constraints, which shows v_L >= 1.
    # The function f(v_L) = v_L - v_L^2 is decreasing for v_L > 0.5.
    # Therefore, to maximize her profit on the feasible interval [1, 15], Chelsea must
    # choose the smallest possible value for v_L.
    v_L = 1.0

    # Step 3: Calculate p_L from v_L
    # Using the binding IR_L constraint: p_L = v_L^2
    p_L = v_L**2
    print(f"Calculation for p_L: p_L = v_L^2 = {v_L:.0f}^2 = {p_L:.0f}")


    # Step 4: Solve for v_H
    # Chelsea's profit does not depend on v_H, so we apply the "no distortion at the top"
    # principle. We choose v_H at its efficient level, which is its maximum possible value
    # as determined by the high type's productivity constraint.
    # Productivity constraint for the high type: v_H <= 30.
    v_H = 30.0

    # Step 5: Calculate p_H from v_H and v_L
    # Using the binding IC_H constraint: p_H = v_H + v_L^2 - v_L
    p_H = v_H + v_L**2 - v_L
    print(f"Calculation for p_H: p_H = v_H + v_L^2 - v_L = {v_H:.0f} + {v_L:.0f}^2 - {v_L:.0f} = {p_H:.0f}")

    print("\nThe separating equilibrium is defined by the following pair of contracts:")
    print(f"Contract for Low Type (v_L, p_L): ({v_L:.0f}, {p_L:.0f})")
    print(f"Contract for High Type (v_H, p_H): ({v_H:.0f}, {p_H:.0f})")

solve_contract_problem()