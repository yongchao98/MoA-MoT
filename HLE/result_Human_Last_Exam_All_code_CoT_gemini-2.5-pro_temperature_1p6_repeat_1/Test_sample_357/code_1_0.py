import sympy as sp

def solve_contract_problem():
    """
    This function solves for the separating equilibrium contracts in the given
    adverse selection problem.
    """
    # Step 1 & 2: Define variables and the Principal's problem
    v_L, p_L, v_H, p_H = sp.symbols('v_L p_L v_H p_H')
    q = sp.Rational(5, 6) # Probability of High type

    # Agent utility functions
    U_H = p_H - v_H
    U_L = p_L - v_L**2

    # Principal's (Chelsea's) expected utility
    E_U_C = q * (v_H - p_H) + (1 - q) * (v_L - p_L)

    print("Step 1: Setting up the problem.")
    print(f"Chelsea's expected utility to maximize: E[U_C] = {E_U_C}\n")

    # Step 3: Assume which constraints bind
    # We hypothesize a non-standard case where the low-type's IC constraint
    # and the high-type's IR constraint bind.
    # (IR-H binds): p_H - v_H = 0
    # (IC-L binds): p_L - v_L^2 = p_H - v_H^2

    print("Step 2: Assuming binding constraints.")
    print("Binding IR-H: p_H - v_H = 0")
    print("Binding IC-L: p_L - v_L^2 = p_H - v_H^2\n")

    # Step 4: Solve for prices
    # From binding IR-H:
    p_H_sol = v_H
    # Substitute into binding IC-L:
    # p_L - v_L^2 = v_H - v_H^2 => p_L = v_L^2 + v_H - v_H^2
    p_L_sol = v_L**2 + v_H - v_H**2

    print("Step 3: Derive prices in terms of values.")
    print(f"From IR-H, we get: p_H = {p_H_sol}")
    print(f"From IC-L, we get: p_L = {p_L_sol}\n")

    # Substitute prices into Chelsea's expected utility
    E_U_C_simplified = E_U_C.subs([(p_H, p_H_sol), (p_L, p_L_sol)])
    E_U_C_simplified = sp.simplify(E_U_C_simplified)
    
    print("Step 4: Formulate Chelsea's simplified optimization problem.")
    print(f"Chelsea's utility becomes: E[U_C] = {E_U_C_simplified}")
    # Let's define g(v) = v - v^2
    g_v_L = v_L - v_L**2
    g_v_H = v_H - v_H**2
    print(f"This can be written as ({1-q}) * (({g_v_L}) - ({g_v_H}))\n")

    # Step 5: Optimize v_L and v_H
    # To maximize E[U_C], Chelsea must:
    # 1. Maximize g(v_L) = v_L - v_L^2
    # 2. Minimize g(v_H) = v_H - v_H^2
    
    # 1. Maximizing g(v_L)
    v_L_opt = sp.solve(sp.diff(g_v_L, v_L), v_L)[0]
    # Check productivity constraint: v_L <= 15. v_L_opt = 0.5 <= 15. This is satisfied.

    # 2. Minimizing g(v_H)
    # We must check the non-binding constraints.
    # (IR-L): p_L - v_L^2 >= 0  => (v_L^2 + v_H - v_H^2) - v_L^2 >= 0 => v_H - v_H^2 >= 0
    # This implies v_H * (1 - v_H) >= 0, so 0 <= v_H <= 1.
    # To minimize g(v_H) on the interval [0, 1], we check the boundaries.
    # g(0) = 0. g(1) = 0. So the minimum is 0, at v_H = 0 or v_H = 1.
    
    print("Step 5: Solving the optimization.")
    print(f"To maximize her utility, Chelsea sets v_L to maximize {g_v_L}.")
    print(f"The maximum occurs at v_L = {v_L_opt}.")
    print(f"She sets v_H to minimize {g_v_H} subject to constraints, which implies v_H must be in [0,1].")
    print("The minimum of g(v_H) on [0,1] occurs at v_H = 0 or v_H = 1.\n")
    
    # Step 6: Select the solution and find the contracts
    # For a separating equilibrium with conventional sorting (v_H > v_L), we choose v_H=1.
    v_H_opt = 1.0

    # Calculate final prices
    p_H_final = v_H_opt
    p_L_final = v_L_opt**2 + v_H_opt - v_H_opt**2

    # Final contract values
    final_v_L = float(v_L_opt)
    final_p_L = float(p_L_final)
    final_v_H = float(v_H_opt)
    final_p_H = float(p_H_final)

    print("Step 6: Determine the final contracts.")
    print("Choosing the conventional sorting v_H > v_L gives v_H = 1.0.")
    print("The final values are:")
    print(f"v_L = {final_v_L}")
    print(f"p_L = {final_p_L}")
    print(f"v_H = {final_v_H}")
    print(f"p_H = {final_p_H}\n")

    print("The final pair of contracts constituting the separating equilibrium is:")
    print(f"((v_L, p_L), (v_H, p_H)) = (({final_v_L}, {final_p_L}), ({final_v_H}, {final_p_H}))")

solve_contract_problem()