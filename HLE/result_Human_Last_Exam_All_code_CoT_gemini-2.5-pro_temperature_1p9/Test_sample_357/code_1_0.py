import sympy

def find_separating_equilibrium():
    """
    This function solves for the separating equilibrium contracts in the given
    adverse selection problem.

    It follows these steps:
    1. Defines all symbolic variables for the model.
    2. Sets up the principal's (Chelsea's) profit function.
    3. Sets up the agents' (handymen's) utility functions.
    4. Assumes the standard binding constraints (IR_L and IC_H) to solve for prices.
    5. Simplifies the principal's profit maximization problem.
    6. Identifies that the IC_L constraint must also bind to maximize profit.
    7. Solves for v_H in terms of v_L for a separating equilibrium.
    8. Substitutes this into the profit function, creating a function of only v_L.
    9. Analyzes this function over the feasible region for v_L (from the IR_H constraint).
    10. Finds the optimal v_L that maximizes profit.
    11. Calculates the final contract values (v_L, p_L, v_H, p_H).
    """
    # 1. Define symbolic variables
    v_L, v_H, p_L, p_H = sympy.symbols('v_L v_H p_L p_H', real=True, nonnegative=True)

    # 2. Define Chelsea's profit function
    pi_H = 30 * v_H - p_H
    pi_L = (20 - sympy.Rational(1, 3) * v_L) * v_L - p_L
    expected_pi = sympy.Rational(5, 6) * pi_H + sympy.Rational(1, 6) * pi_L

    # 3. Define handyman utilities
    u_H_vL_pL = p_L - v_L
    u_H_vH_pH = p_H - v_H
    u_L_vL_pL = p_L - v_L**2
    u_L_vH_pH = p_H - v_H**2

    # 4. Assume IR_L and IC_H constraints bind and solve for prices
    p_L_sol = sympy.solve(sympy.Eq(u_L_vL_pL, 0), p_L)[0]
    p_H_sol = sympy.solve(sympy.Eq(u_H_vH_pH, u_H_vL_pL).subs(p_L, p_L_sol), p_H)[0]

    # 5. Substitute prices into profit function
    profit_in_v = expected_pi.subs([(p_L, p_L_sol), (p_H, p_H_sol)])
    
    # 6. From profit maximization, the IC_L constraint must also bind
    ic_l_constraint = u_L_vL_pL.subs(p_L, p_L_sol) >= u_L_vH_pH.subs(p_H, p_H_sol)
    ic_l_binding = sympy.Eq(ic_l_constraint.lhs, ic_l_constraint.rhs)
    
    # 7. Solve for v_H. There are two solutions: v_H=v_L (pooling) and v_H=1-v_L (separating)
    v_H_expr = 1 - v_L

    # 8. Substitute v_H into the profit function
    profit_in_vL = profit_in_v.subs(v_H, v_H_expr)

    # 9. Analyze the profit function over the feasible set for v_L
    # The IR_H constraint (p_H - v_H >= 0) implies v_L*(v_L-1) >= 0.
    # The feasible set for v_L is {0} U [1, infinity).
    # The derivative of the profit function w.r.t v_L is negative for v_L >= 0,
    # meaning profit is a decreasing function of v_L.
    
    # 10. To maximize a decreasing function, choose the minimum value in the feasible set.
    v_L_final = 0

    # 11. Calculate the final contract values
    v_H_final = v_H_expr.subs(v_L, v_L_final)
    p_L_final = p_L_sol.subs(v_L, v_L_final)
    p_H_final = p_H_sol.subs([(v_L, v_L_final), (v_H, v_H_final)])
    
    # The final equations for the contracts are (v_L, p_L) and (v_H, p_H)
    # I will output each number for these final contract pairs.
    print("The separating equilibrium constitutes the following pair of contracts:")
    print(f"Contract offered to Low-Type Handyman: (v_L, p_L) = ({v_L_final}, {p_L_final})")
    print(f"Contract offered to High-Type Handyman: (v_H, p_H) = ({v_H_final}, {p_H_final})")


find_separating_equilibrium()