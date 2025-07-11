import math

def solve_for_contracts():
    """
    This function calculates the contracts (v_L, p_L) and (v_H, p_H) for the separating equilibrium.
    """
    # Step 1: Solve the binding IR_L constraint for v_L.
    # The equation is v_L^2 + (1/3)v_L - 20 = 0.
    # Using the quadratic formula: v_L = (-b +/- sqrt(b^2 - 4ac)) / 2a
    # a=1, b=1/3, c=-20
    a = 1
    b = 1/3
    c = -20
    
    # Calculate the discriminant
    discriminant = b**2 - 4*a*c
    
    # v_L must be non-negative (as an effort level)
    v_L = (-b + math.sqrt(discriminant)) / (2*a)
    
    # Step 2: Calculate p_L using the zero-profit condition.
    # p_L = 20 - (1/3)*v_L.
    # Since IR_L is binding, p_L is also equal to v_L^2.
    p_L = v_L**2
    
    # Step 3: Solve the binding IC_H constraint for v_H.
    # 30 - v_H = p_L - v_L
    # v_H = 30 - p_L + v_L
    v_H = 30 - p_L + v_L
    
    # Step 4: The price for the high type is fixed by the zero-profit condition.
    p_H = 30.0

    # Print the resulting contract specifications
    print("The separating equilibrium is defined by the following pair of contracts:")
    print("\nContract for Low Type (L):")
    print(f"(v_L, p_L) = ({v_L}, {p_L})")
    print("\nContract for High Type (H):")
    print(f"(v_H, p_H) = ({v_H}, {p_H})")

    # Final check of the binding IC_H constraint to show consistency
    # u_H(C_H) = p_H - v_H
    # u_H(C_L) = p_L - v_L
    # print(f"\nCheck: u_H at H-contract: {p_H - v_H}, u_H at L-contract: {p_L - v_L}")

solve_for_contracts()