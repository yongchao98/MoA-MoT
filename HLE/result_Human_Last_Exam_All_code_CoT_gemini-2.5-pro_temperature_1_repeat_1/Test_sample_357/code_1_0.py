import sympy

def solve_contracts():
    """
    This function calculates the separating equilibrium contracts based on the problem's parameters.
    """
    # Define variables using sympy for symbolic manipulation if needed
    v_L_sym = sympy.symbols('v_L')

    # Define productivity and price functions based on the zero-profit condition
    theta_H = 30
    def theta_L_func(v):
        return 20 - (1/3) * v

    def p_H_func(v):
        return theta_H * v

    def p_L_func(v):
        return theta_L_func(v) * v

    # Step 1: Find the value v_L for the Low-type contract.
    # We maximize the High-type's utility on the Low-type's contract, subject to
    # the Low-type's Individual Rationality (IR) constraint.
    # IR-L: u_L = p_L - v_L^2 >= 0 => 20*v_L - (4/3)*v_L^2 >= 0 => 0 <= v_L <= 15
    # The range for v_L is [0, 15].
    # Utility for H-type if they take L-type's contract: u_H_on_L = p_L(v_L) - v_L
    u_H_on_L_expr = p_L_func(v_L_sym) - v_L_sym

    # To maximize this on [0, 15], we check the derivative.
    # deriv = 19 - (2/3)*v_L, which is positive on [0, 15].
    # Thus, the maximum is at the upper bound.
    v_L_sol = 15.0

    # Calculate p_L for this v_L
    p_L_sol = p_L_func(v_L_sol)

    # Step 2: Find the value v_H for the High-type contract.
    # This requires satisfying two incentive compatibility (IC) constraints.
    # IC-H constraint: u_H(v_H) >= u_H(v_L_sol)
    # 29*v_H >= p_L_sol - v_L_sol
    u_H_from_L_contract = p_L_sol - v_L_sol
    min_v_H_from_ICH = u_H_from_L_contract / 29

    # IC-L constraint: u_L(v_L_sol) >= u_L(v_H)
    # 0 >= p_H(v_H) - v_H^2  => v_H * (v_H - 30) >= 0
    # This implies v_H <= 0 or v_H >= 30.

    # Combine constraints: v_H >= min_v_H_from_ICH (approx 7.24) AND (v_H <= 0 or v_H >= 30)
    # This leads to the condition v_H >= 30.
    # The equilibrium forms at the lowest possible value for v_H that satisfies this.
    v_H_sol = 30.0

    # Calculate p_H for this v_H
    p_H_sol = p_H_func(v_H_sol)
    
    # Output the results
    print("The separating equilibrium consists of the following pair of contracts:")
    print("\nContract for Low-Type Handyman (v_L, p_L):")
    print(f"Value added (v_L) = {v_L_sol}")
    print(f"Price (p_L) = {p_L_sol}")
    
    print("\nContract for High-Type Handyman (v_H, p_H):")
    print(f"Value added (v_H) = {v_H_sol}")
    print(f"Price (p_H) = {p_H_sol}")

solve_contracts()