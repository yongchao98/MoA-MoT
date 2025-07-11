import sympy

def find_separating_equilibrium():
    """
    Calculates the pair of contracts for a separating equilibrium.
    """
    # Step 1: Define constraints and utilities based on the problem statement.
    
    # Feasibility constraint for High Type
    v_h_max = 30
    
    # Feasibility constraint for Low Type: v <= 20 - (1/3)*v => (4/3)*v <= 20
    v_l_max = 20 / (4/3)

    # Step 2: Determine the Low Type's optimal contract.
    # Under competition, profit is zero, so p = v.
    # The Low Type chooses v to maximize utility u_L = v - v^2, subject to v <= v_l_max.
    v = sympy.Symbol('v')
    u_l = v - v**2
    
    # Find the maximum by taking the derivative and solving for v.
    du_l_dv = sympy.diff(u_l, v)
    # Solve 1 - 2v = 0
    v_l_optimal_unconstrained = sympy.solve(du_l_dv, v)[0]

    # Check if the unconstrained optimum is within the feasible range.
    if v_l_optimal_unconstrained <= v_l_max:
        v_l = v_l_optimal_unconstrained
    else:
        # If the optimum is not feasible, the choice will be at the boundary.
        # (This is not the case here).
        v_l = v_l_max
    
    # The price p_l is equal to the value v_l due to competition.
    p_l = v_l

    # Step 3: Determine the High Type's optimal contract.
    # The High Type chooses v to maximize utility u_H = v - v.
    # The utility is always 0, so the agent is indifferent.
    # We apply the "no distortion at the top" principle: select the max feasible value.
    v_h = v_h_max
    
    # The price p_h is equal to the value v_h.
    p_h = v_h

    # Step 4: Output the results.
    print("A separating equilibrium is constituted by the following pair of contracts:")
    print("Contract for Low Type (v_L, p_L):")
    # Using .evalf() to ensure floating point representation for clarity
    print(f"(v_L = {v_l.evalf()}, p_L = {p_l.evalf()})")
    
    print("\nContract for High Type (v_H, p_H):")
    print(f"(v_H = {v_h}, p_H = {p_h})")

# Run the function to find and print the answer.
find_separating_equilibrium()