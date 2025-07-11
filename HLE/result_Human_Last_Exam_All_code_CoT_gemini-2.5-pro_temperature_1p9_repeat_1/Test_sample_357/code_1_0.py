import sympy

def solve_handyman_contracts():
    """
    Calculates the separating equilibrium contracts for Chelsea and the handymen.
    """

    # Step 1 & 2: Define handyman utilities under the zero-profit condition (p=v)
    # due to competition among homeowners.
    
    # Low-type's utility to be maximized: u_L = v - v^2
    v = sympy.Symbol('v')
    u_L = v - v**2
    
    # Step 3: Determine the low-type's optimal contract.
    # Maximize u_L by finding where its derivative is zero.
    du_L_dv = sympy.diff(u_L, v)
    
    # Solve d(u_L)/dv = 0 for v
    v_L_solution = sympy.solve(du_L_dv, v)[0]
    v_L = float(v_L_solution)

    # Check the low-type's productivity constraint: v_L <= 20 - (1/3)*v_L
    # This simplifies to v_L <= 15
    productivity_limit_L = 15
    if v_L > productivity_limit_L:
        # This case is not reached, but good practice to handle.
        v_L = productivity_limit_L

    # The price p_L is equal to the value v_L.
    p_L = v_L
    
    # Step 4: Determine the high-type's contract.
    # The high-type is indifferent and chooses their maximum productive value.
    theta_H = 30
    v_H = theta_H
    
    # The price p_H is equal to the value v_H.
    p_H = v_H
    
    # Print the results and the logic
    print("Based on the analysis, the separating equilibrium contracts are:")
    print("-" * 50)
    print("Contract for the Low-Type Handyman (L):")
    print(f"The utility function to maximize is u_L(v) = v - v^2.")
    print(f"The optimal value is found by solving d(u_L)/dv = {du_L_dv} = 0, which gives v_L = {v_L}")
    print(f"The resulting contract (v_L, p_L) is ({v_L}, {p_L}).")
    print("-" * 50)
    print("Contract for the High-Type Handyman (H):")
    print(f"The high-type's utility u_H(v) = v - v is always 0, so they are indifferent.")
    print(f"They choose their maximum productive capacity, v_H = {v_H}.")
    print(f"The resulting contract (v_H, p_H) is ({v_H}, {p_H}).")
    print("-" * 50)

    # Final Answer Formatting
    print("The final pair of contracts constituting the separating equilibrium is:")
    print(f"((v_L, p_L), (v_H, p_H)) = (({v_L}, {p_L}), ({v_H}, {p_H}))")

solve_handyman_contracts()