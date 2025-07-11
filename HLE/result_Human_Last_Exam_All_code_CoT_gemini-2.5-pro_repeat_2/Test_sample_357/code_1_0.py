import sympy

def solve_contract_problem():
    """
    This script solves for the separating equilibrium contracts in the given principal-agent problem.
    It follows the standard method of assuming which constraints are binding and then maximizing
    the principal's utility function.
    """
    # Define symbolic variables for the contracts
    v_L, p_L, v_H, p_H = sympy.symbols('v_L p_L v_H p_H')
    
    print("Step 1: Define the constraints for a separating equilibrium.")
    print("IR-L: p_L - v_L**2 >= 0")
    print("IR-H: p_H - v_H >= 0")
    print("IC-L: p_L - v_L**2 >= p_H - v_H**2")
    print("IC-H: p_H - v_H >= p_L - v_L")
    print("-" * 30)

    print("Step 2: Assume IR-L and IC-H are binding to express prices in terms of values.")
    # From IR-L binding: p_L - v_L^2 = 0
    price_L_expr = v_L**2
    print(f"From IR-L binding, we get: p_L = {price_L_expr}")

    # From IC-H binding: p_H - v_H = p_L - v_L
    price_H_expr = v_H + p_L - v_L
    price_H_expr = price_H_expr.subs(p_L, price_L_expr)
    print(f"From IC-H binding, we get: p_H = v_H + p_L - v_L. Substituting for p_L gives: p_H = {price_H_expr}")
    print("-" * 30)

    print("Step 3: Formulate Chelsea's expected profit function.")
    # Probabilities and productivities
    prob_H = sympy.Rational(5, 6)
    prob_L = sympy.Rational(1, 6)
    theta_H = 30
    theta_L = 20 - sympy.Rational(1, 3) * v_L

    profit = prob_H * (theta_H - p_H) + prob_L * (theta_L - p_L)
    profit = profit.subs(p_L, price_L_expr).subs(p_H, price_H_expr)
    
    print(f"Chelsea's Expected Profit E[pi] = {prob_H}*(30 - p_H) + {prob_L}*((20 - (1/3)*v_L) - p_L)")
    print("Substituting prices, the profit function to maximize is:")
    simplified_profit = sympy.simplify(profit)
    print(f"E[pi] = {simplified_profit}")
    print("-" * 30)
    
    print("Step 4: Check the remaining constraints (IR-H and IC-L) to find the feasible range for v_L and v_H.")
    # IR-H constraint: p_H - v_H >= 0
    ir_h_constraint = price_H_expr - v_H
    print(f"IR-H becomes: {ir_h_constraint} >= 0, which simplifies to v_L*(v_L - 1) >= 0.")
    print("Assuming v_L >= 0, this means v_L must be 0 or v_L >= 1.")
    
    # IC-L constraint: p_L - v_L^2 >= p_H - v_H^2
    ic_l_constraint_lhs = price_L_expr - v_L**2
    ic_l_constraint_rhs = price_H_expr - v_H**2
    print(f"IC-L becomes: {ic_l_constraint_lhs} >= {ic_l_constraint_rhs}, which simplifies to v_H**2 - v_H >= v_L**2 - v_L.")
    print("-" * 30)

    print("Step 5: Maximize the profit function by choosing v_L and v_H.")
    # The profit function is E[pi] = 85/3 - 5*v_H/6 - v_L**2 + 7*v_L/9
    # To maximize, we need to minimize v_H and maximize the v_L part: g(v_L) = -v_L**2 + 7*v_L/9.
    
    # Analyze g(v_L) over the feasible range for v_L (v_L=0 or v_L>=1)
    g_vL = -v_L**2 + sympy.Rational(7, 9) * v_L
    g_at_0 = g_vL.subs(v_L, 0)
    g_at_1 = g_vL.subs(v_L, 1)
    
    print(f"The part of the profit depending on v_L is g(v_L) = {g_vL}.")
    print(f"At v_L=0, g(0) = {g_at_0}.")
    print(f"At v_L=1, g(1) = {g_at_1}.")
    print("For v_L > 1, g(v_L) is decreasing. Thus, g(v_L) is maximized at v_L=0.")
    
    optimal_v_L = 0
    print(f"Optimal v_L = {optimal_v_L}")

    # With v_L=0, find optimal v_H.
    # Profit becomes E[pi] = 85/3 - 5*v_H/6. To maximize, minimize v_H.
    # The constraint on v_H is v_H**2 - v_H >= 0, which means v_H>=1 (for a separating equilibrium where v_H!=v_L).
    optimal_v_H = 1
    print(f"With v_L=0, the IC-L constraint on v_H is v_H*(v_H-1) >= 0. For a separating equilibrium, v_H must be different from v_L=0, so v_H >= 1.")
    print(f"To maximize profit, we must minimize v_H. So, optimal v_H = {optimal_v_H}.")
    print("-" * 30)
    
    print("Step 6: Calculate the final contracts based on optimal values.")
    final_v_L = 0
    final_p_L_val = final_v_L**2
    
    final_v_H = 1
    final_p_H_val = final_v_H + final_v_L**2 - final_v_L
    
    print("The separating equilibrium contract for the low type (L) is:")
    print(f"(v_L, p_L) = ({final_v_L}, {final_p_L_val})")
    print(f"The price p_L is calculated as: p_L = v_L^2 = {final_v_L}^2 = {final_p_L_val}")
    
    print("\nThe separating equilibrium contract for the high type (H) is:")
    print(f"(v_H, p_H) = ({final_v_H}, {final_p_H_val})")
    print(f"The price p_H is calculated as: p_H = v_H + v_L^2 - v_L = {final_v_H} + {final_v_L}^2 - {final_v_L} = {final_p_H_val}")

solve_contract_problem()