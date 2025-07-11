import sympy
from sympy import symbols, solve, diff, Max

def solve_contract_problem():
    """
    Solves for the separating equilibrium contracts in the given problem.
    """
    # Define symbolic variables
    v_H, v_L, p_H, p_L = symbols('v_H v_L p_H p_L', real=True)
    q = sympy.Rational(5, 6) # Probability of High type

    # From the problem's logic (standard screening model), we assume:
    # 1. Low type's Individual Rationality (IR_L) constraint is binding.
    #    p_L - v_L**2 = 0  => p_L = v_L**2
    p_L_sol = v_L**2

    # 2. High type's Incentive Compatibility (IC_H) constraint is binding.
    #    p_H - v_H = p_L - v_L
    #    Substituting p_L gives the information rent for the high type.
    #    p_H = v_H + p_L - v_L
    p_H_sol = v_H + p_L_sol - v_L

    # Now, check the remaining two constraints (IR_H and IC_L)
    # IR_H check: p_H - v_H >= 0 => (v_H + v_L**2 - v_L) - v_H >= 0 => v_L**2 - v_L >= 0
    # This implies v_L <= 0 or v_L >= 1. Assuming v >= 0, this means v_L=0 or v_L>=1.

    # IC_L check: p_L - v_L**2 >= p_H - v_H**2
    # Since IR_L is binding, p_L - v_L**2 = 0.
    # So, 0 >= p_H - v_H**2 => 0 >= (v_H + v_L**2 - v_L) - v_H**2
    # This simplifies to: v_H**2 - v_H <= v_L**2 - v_L

    # Chelsea's expected profit function
    profit_H = v_H - p_H_sol
    profit_L = v_L - p_L_sol
    expected_profit = q * profit_H + (1 - q) * profit_L

    # Substitute the simplified profits.
    # The term v_H cancels out in the high-type profit calculation, which is a known issue.
    # Let's re-evaluate Chelsea's objective function based on agent utilities (rents).
    # u_L = 0 (from binding IR_L)
    # u_H = p_H - v_H = v_L**2 - v_L (from binding IC_H). This is the info rent.
    # Profit from H: v_H - p_H = v_H - (v_H + u_H) = -u_H = -(v_L**2 - v_L) = v_L - v_L**2
    # Profit from L: v_L - p_L = v_L - v_L**2
    # This implies profit is independent of type once contracts are set. This is unusual.
    
    # A more careful formulation:
    # Let u_L = p_L - v_L**2 and u_H = p_H - v_H.
    # Chelsea's problem: max E[pi] = q*(v_H - (u_H+v_H)) + (1-q)*(v_L - (u_L+v_L**2))
    # E[pi] = -q*u_H + (1-q)*(v_L - u_L - v_L**2)
    # Chelsea wants to minimize u_H and u_L.
    
    # Minimal rents subject to all constraints:
    u_L_min = 0 # From IR_L
    u_H_min = Max(0, v_L**2 - v_L) # From IR_H and IC_H
    
    # Substitute minimal rents into profit function
    E_pi_expr = -q * Max(0, v_L**2 - v_L) + (1-q) * (v_L - 0 - v_L**2)
    
    # We maximize this profit w.r.t v_L and v_H, subject to the constraints.
    # Constraint: v_H**2 - v_H >= u_H_min + v_L**2 - v_L (derived from IC_L after substitution, incorrect, re-deriving)
    # Correct IC_L check: u_L >= u_H + v_H - v_H**2 => 0 >= u_H_min + v_H - v_H**2 => u_H_min <= v_H**2 - v_H
    
    # Case 1: v_L is in [0, 1]. Then v_L**2 - v_L <= 0, so Max(0, v_L**2 - v_L) = 0.
    # E_pi_case1 = (1-q)*(v_L - v_L**2). To maximize, we take the derivative w.r.t v_L.
    f_vL = v_L - v_L**2
    df_vL = diff(f_vL, v_L) # 1 - 2*v_L
    v_L_opt_case1 = solve(df_vL, v_L)[0] # Gives 1/2
    
    # At v_L=1/2, the condition v_L in [0,1] is met.
    # The constraint on v_H is u_H_min <= v_H**2 - v_H => 0 <= v_H**2 - v_H => v_H(v_H-1)>=0.
    # This requires v_H <= 0 or v_H >= 1. Since v is value, v_H=0 or v_H>=1.
    # Productivity constraints: v_L <= 15 and v_H <= 30 are satisfied.
    # For separation, we need (v_H,p_H) != (v_L,p_L).
    
    # Case 2: v_L > 1. Then v_L**2 - v_L > 0.
    # E_pi_case2 = -q*(v_L**2 - v_L) + (1-q)*(v_L-v_L**2) = v_L - v_L**2.
    # For v_L > 1, this function is negative and decreasing.
    # Max profit in this region is approached as v_L -> 1, where profit approaches 0.
    
    # Comparing cases, the max profit is in Case 1.
    v_L_final = v_L_opt_case1

    # We can choose any v_H satisfying v_H >= 1 and v_H <= 30. A canonical choice is the lowest possible value.
    v_H_final = 1
    
    # Now we find the prices using the binding constraint formulas
    p_L_final = v_L_final**2
    p_H_final = v_H_final + v_L_final**2 - v_L_final

    print("The separating equilibrium consists of the following two contracts (v, p):")
    print("\nContract for the Low Type (L):")
    print(f"Value v_L = {v_L_final}")
    print(f"Price p_L = {p_L_final}")
    print(f"(v_L, p_L) = ({v_L_final}, {p_L_final})")
    
    print("\nContract for the High Type (H):")
    print(f"Value v_H = {v_H_final}")
    print(f"Price p_H = {p_H_final}")
    print(f"(v_H, p_H) = ({v_H_final}, {p_H_final})")

solve_contract_problem()