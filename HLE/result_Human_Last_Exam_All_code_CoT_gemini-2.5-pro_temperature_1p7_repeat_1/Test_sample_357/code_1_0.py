import sympy

def solve_contract_problem():
    """
    This function solves for the separating equilibrium contracts in the
    handyman screening problem.
    """
    # Step 1: Define symbols for our variables
    v_L, p_L, v_H, p_H = sympy.symbols('v_L p_L v_H p_H')

    # Step 2: Set up the constraints and objective function based on the plan.
    # From economic theory, two constraints will be binding:
    # 1. The Low type's Individual Rationality (IR_L) constraint.
    #    p_L - v_L**2 = 0
    # 2. The High type's Incentive Compatibility (IC_H) constraint.
    #    p_H - v_H = p_L - v_L

    # From these binding constraints, we can express prices in terms of values.
    eq_p_L = sympy.Eq(p_L, v_L**2)
    p_L_sol = sympy.solve(eq_p_L, p_L)[0]

    eq_p_H = sympy.Eq(p_H - v_H, p_L_sol - v_L)
    p_H_sol = sympy.solve(eq_p_H, p_H)[0]
    
    # Step 3: Analyze the remaining (non-binding) constraints to find bounds on v_L and v_H.
    # - IR_H (p_H - v_H >= 0) implies v_L**2 - v_L >= 0. Since v_L must be positive, this means v_L >= 1.
    # - IC_L (p_L - v_L**2 >= p_H - v_H**2) implies v_H**2 - v_H >= v_L**2 - v_L. Since f(x)=x**2-x is increasing for x>0.5, this means v_H >= v_L.
    # - Productivity for L-type (v_L <= 20 - (1/3)*v_L) implies (4/3)*v_L <= 20, so v_L <= 15.
    # - Productivity for H-type (v_H <= 30).

    # Step 4: Maximize Chelsea's expected profit function.
    # E[Profit] = (5/6)*(v_H - p_H) + (1/6)*(v_L - p_L)
    # Substituting p_H and p_L from the binding constraints:
    # E[Profit] = (5/6)*(v_H - (v_H + v_L**2 - v_L)) + (1/6)*(v_L - v_L**2)
    # This simplifies to: E[Profit] = v_L - v_L**2

    # We need to maximize f(v_L) = v_L - v_L**2 subject to 1 <= v_L <= 15.
    # The function's maximum is at v_L = 0.5. On the interval [1, 15], the function is decreasing.
    # Therefore, Chelsea maximizes her profit by choosing the smallest possible v_L.
    v_L_optimal = 1.0

    # Step 5: Determine v_H.
    # Chelsea's profit doesn't depend on v_H. We apply the "no distortion at the top"
    # principle, setting v_H to its maximum feasible value to maximize production.
    # The constraints on v_H are v_H >= v_L and v_H <= 30.
    # v_H_optimal must be in [1.0, 30.0]. We pick the maximum value.
    v_H_optimal = 30.0

    # For this to be a separating equilibrium, the contracts must differ.
    # (v_L, p_L) will be (1, 1). (v_H, p_H) will be (30, 30), which is different.

    # Step 6: Calculate the final prices for the contracts.
    p_L_optimal = v_L_optimal**2
    p_H_optimal = v_H_optimal + p_L_optimal - v_L_optimal

    # Step 7: Print the final contract parameters.
    print("The separating equilibrium is defined by two contracts: one for the Low type (L) and one for the High type (H).")
    print("\nContract for Low Type Handyman (L):")
    print(f"Value added (v_L) = {v_L_optimal}")
    print(f"Price paid (p_L) = {p_L_optimal}")

    print("\nContract for High Type Handyman (H):")
    print(f"Value added (v_H) = {v_H_optimal}")
    print(f"Price paid (p_H) = {p_H_optimal}")

solve_contract_problem()