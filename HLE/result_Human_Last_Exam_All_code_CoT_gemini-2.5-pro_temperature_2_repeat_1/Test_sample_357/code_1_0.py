import math

def solve_separating_equilibrium():
    """
    This function calculates the separating equilibrium contracts for Chelsea
    by following a step-by-step derivation and verification process.
    """
    pi_H = 5/6
    pi_L = 1 - pi_H

    # Step 1: Theoretical Derivations (as outlined in the plan)
    # Based on the analysis, we determined that in this specific problem,
    # the IR-H and IC-L constraints are binding.
    # IR-H binds: p_H - v_H = 0  => p_H = v_H
    # IC-L binds: p_L - v_L^2 = p_H - v_H^2
    #
    # Substituting p_H = v_H into the second binding constraint gives:
    # p_L - v_L^2 = v_H - v_H^2 => p_L = v_L^2 + v_H - v_H^2
    #
    # The principal's expected utility to maximize is:
    # E_U = pi_L * (v_L - p_L) + pi_H * (v_H - p_H)
    # E_U = pi_L * (v_L - (v_L^2 + v_H - v_H^2)) + pi_H * (v_H - v_H)
    # E_U = pi_L * (v_L - v_L^2 - v_H + v_H^2)
    #
    # We must maximize this utility subject to the derived constraints from
    # the non-binding IR and IC conditions.
    # From non-binding IR-L: v_H - v_H^2 >= 0 => 0 <= v_H <= 1
    # From non-binding IC-H: v_L - v_L^2 >= v_H - v_H^2.
    # Let h(v) = v - v^2. For v > 0.5, h(v) is decreasing, so h(v_L) >= h(v_H) implies v_L <= v_H.
    # The optimization leads to a corner solution. Maximizing (v_L - v_L^2) and minimizing
    # (v_H - v_H^2) subject to v_L < v_H and 0 <= v_H <= 1 points to:
    # v_L is where v-v^2 is max (v=0.5)
    # v_H is where v-v^2 is min on [0,1] (v=0 or v=1).
    # To maintain v_L < v_H, we must choose v_H=1.

    # Step 2: Calculate the optimal contract parameters based on the derivation
    v_L = 0.5
    v_H = 1.0

    # Step 3: Calculate the prices based on the binding constraints
    p_H = v_H
    p_L = v_L**2 + v_H - v_H**2

    print("--- Separating Equilibrium Contracts ---")
    print(f"Contract for Low Type (v_L, p_L):")
    print(f"  Value v_L = {v_L}")
    print(f"  Price p_L = {p_L}\n")

    print(f"Contract for High Type (v_H, p_H):")
    print(f"  Value v_H = {v_H}")
    print(f"  Price p_H = {p_H}\n")

    # Step 4: Verify the solution by checking all constraints
    print("--- Verification of Constraints ---")

    # Individual Rationality (IR)
    ir_L = p_L - v_L**2
    ir_H = p_H - v_H
    print(f"(IR-L) Low Type's Utility: p_L - v_L^2 = {p_L} - {v_L}^2 = {ir_L:.4f} >= 0.  ({ir_L >= -1e-9})")
    print(f"(IR-H) High Type's Utility: p_H - v_H = {p_H} - {v_H} = {ir_H:.4f} >= 0.  ({ir_H >= -1e-9})\n")

    # Incentive Compatibility (IC)
    # Low type should prefer (v_L, p_L)
    utility_L_L = p_L - v_L**2
    utility_L_H = p_H - v_H**2
    ic_L_check = utility_L_L >= utility_L_H
    print("(IC-L) Low type prefers (v_L, p_L):")
    print(f"  Utility from (v_L,p_L): {utility_L_L:.4f}")
    print(f"  Utility from (v_H,p_H): {utility_L_H:.4f}")
    print(f"  Is {utility_L_L:.4f} >= {utility_L_H:.4f}? {ic_L_check}\n")

    # High type should prefer (v_H, p_H)
    utility_H_H = p_H - v_H
    utility_H_L = p_L - v_L
    ic_H_check = utility_H_H >= utility_H_L
    print("(IC-H) High type prefers (v_H, p_H):")
    print(f"  Utility from (v_H,p_H): {utility_H_H:.4f}")
    print(f"  Utility from (v_L,p_L): {utility_H_L:.4f}")
    print(f"  Is {utility_H_H:.4f} >= {utility_H_L:.4f}? {ic_H_check}\n")

    # Feasibility Constraints from productivity
    theta_H = 30
    theta_L = 20 - (1/3) * v_L
    feasibility_H = v_H <= theta_H
    feasibility_L = v_L <= theta_L
    print("--- Verification of Feasibility Constraints ---")
    print(f"High type: v_H <= 30?  {v_H} <= 30 is {feasibility_H}")
    print(f"Low type: v_L <= 20 - v_L/3?  {v_L} <= {theta_L:.4f} is {feasibility_L}\n")

    # Final Answer Formatting
    final_answer = f"<<<(v_L, p_L) = ({v_L}, {p_L}), (v_H, p_H) = ({v_H}, {p_H})>>>"
    print(final_answer)


solve_separating_equilibrium()