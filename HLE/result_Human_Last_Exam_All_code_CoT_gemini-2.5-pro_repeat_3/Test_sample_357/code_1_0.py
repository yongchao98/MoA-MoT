import sympy

def find_separating_equilibrium():
    """
    Solves for the separating equilibrium contracts in the given principal-agent problem.

    This function follows a standard procedure for solving adverse selection problems:
    1. Set up the principal's (Chelsea's) profit maximization problem with IR and IC constraints.
    2. Assume the low-type's IR constraint and the high-type's IC constraint are binding.
    3. Solve for prices (p_L, p_H) in terms of value-added (v_L, v_H).
    4. Substitute prices into the profit function and simplify.
    5. Maximize the resulting profit function with respect to v_L and v_H subject to the remaining constraints.
    6. Calculate the specific contract values.
    """
    # Step 1: Define symbolic variables
    v_H, p_H, v_L, p_L = sympy.symbols('v_H p_H v_L p_L')
    pi_H = sympy.Rational(5, 6)
    pi_L = sympy.Rational(1, 6)

    # Chelsea's expected utility (profit)
    chelsea_profit = pi_H * (v_H - p_H) + pi_L * (v_L - p_L)

    print("Step 1: Define the Principal's objective function.")
    print(f"Chelsea's Expected Profit E[Π] = {pi_H}*(v_H - p_H) + {pi_L}*(v_L - p_L)\n")

    # Step 2 & 3: Assume binding constraints and solve for prices
    # The IR constraint for the low type (L) is u_L >= 0 => p_L - v_L**2 >= 0.
    # We assume it's binding: p_L = v_L**2
    p_L_sol = v_L**2
    print("Step 2: Assume the Low-Type's Individual Rationality (IR-L) constraint is binding.")
    print(f"From p_L - v_L^2 = 0, we get: p_L = {p_L_sol}\n")

    # The IC constraint for the high type (H) is u_H(H) >= u_H(L) => p_H - v_H >= p_L - v_L.
    # We assume it's binding: p_H - v_H = p_L - v_L
    p_H_sol = v_H + p_L_sol - v_L
    print("Step 3: Assume the High-Type's Incentive Compatibility (IC-H) constraint is binding.")
    print(f"From p_H - v_H = p_L - v_L, we get: p_H = v_H + p_L - v_L = {p_H_sol}\n")

    # Step 4: Substitute prices into the profit function
    profit_in_v = chelsea_profit.subs([(p_L, p_L_sol), (p_H, p_H_sol)])
    simplified_profit = sympy.simplify(profit_in_v)
    print("Step 4: Substitute prices into the profit function to express it in terms of v_L and v_H.")
    print(f"E[Π] = {profit_in_v}")
    print(f"Simplified E[Π] = {simplified_profit}\n")

    # Step 5: Optimize with respect to v_L and v_H
    # The profit function E[Π] = v_L - v_L^2 depends only on v_L.
    # We must check the remaining constraints to find the feasible range for v_L.
    # Constraint 1 (IR-H): p_H - v_H >= 0 => (v_H + v_L**2 - v_L) - v_H >= 0 => v_L**2 - v_L >= 0.
    # This implies v_L(v_L - 1) >= 0. For a meaningful repair, v_L > 0, so v_L >= 1.
    # Constraint 2 (Feasibility-L): v_L <= 20 - (1/3)*v_L => (4/3)*v_L <= 20 => v_L <= 15.
    # To maximize f(v_L) = v_L - v_L^2 on the interval v_L in [1, 15], we note that the function
    # is decreasing for v_L > 0.5. Thus, the maximum is at the lowest possible value, v_L = 1.
    v_L_val = 1
    print("Step 5: Optimize the profit function.")
    print("The profit E[Π] = v_L - v_L^2 must be maximized subject to the remaining constraints.")
    print("The IR-H and Feasibility constraints imply that 1 <= v_L <= 15.")
    print("To maximize v_L - v_L^2 in this range, we must choose the smallest value, so v_L = 1.\n")

    # The choice of v_H does not affect profit. This is the "no distortion at the top" principle.
    # The efficient choice of v_H is the one that maximizes surplus, but here surplus is v_H - v_H = 0 for any v_H.
    # The standard convention is to set v_H to its maximum feasible level to resolve this.
    # Feasibility-H: v_H <= 30.
    v_H_val = 30
    print("The choice of v_H is determined by the 'no distortion at the top' principle.")
    print("We set v_H to its maximum feasible value, v_H = 30.\n")

    # Step 6: Calculate final contract parameters
    p_L_val = v_L_val**2
    p_H_val = v_H_val + v_L_val**2 - v_L_val

    print("Step 6: Calculate the final contract parameters based on optimal v_L and v_H.")
    print("The separating equilibrium is the pair of contracts ((v_L, p_L), (v_H, p_H)):")
    print(f"Low-type contract (v_L, p_L):")
    print(f"v_L = {v_L_val}")
    print(f"p_L = v_L^2 = {v_L_val}^2 = {p_L_val}\n")
    print(f"High-type contract (v_H, p_H):")
    print(f"v_H = {v_H_val}")
    print(f"p_H = v_H + v_L^2 - v_L = {v_H_val} + {v_L_val}^2 - {v_L_val} = {p_H_val}\n")
    
    print("The final pair of contracts is ((1, 1), (30, 30)).")

find_separating_equilibrium()