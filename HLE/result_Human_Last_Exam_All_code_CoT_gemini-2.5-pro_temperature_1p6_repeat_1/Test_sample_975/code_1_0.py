import sympy as sp

def display_magnetic_field_solution():
    """
    This function displays the symbolic solution for the magnetic field H
    in the two regions of the specified magnetostatics problem.
    The derivation is based on solving Laplace's equation for the magnetic
    scalar potential with the appropriate boundary conditions.
    """

    # Define symbolic variables
    M0, R, Rp, r, theta = sp.symbols('M_0 R R_p r theta')
    i_r = sp.Symbol('î_r')
    i_theta = sp.Symbol('î_θ')
    
    print("Based on the derivation, the correct option is B.\n")
    print("The magnetic field H(r, theta) is determined as follows:")
    print("-" * 50)

    # --- Region 1: Inside the shield (0 < r < R_p) ---
    print("In the region 0 < r < R_p:")
    
    # Coefficient for the field inside the shield
    H1_coeff_val = M0 * (2*Rp**3 + R**3) / (3*R**3)
    
    # H-field vector components inside the shield
    H1_r_comp_val = -sp.cos(theta)
    H1_theta_comp_val = sp.sin(theta)

    # We will print the terms in the final equation.
    # Printing the scalar coefficient part
    print(f"H = (M_0 * (2*R_p**3 + R**3) / (3*R**3)) * ( ... )")
    
    # Constructing the vector part as a string to show the structure
    H1_vector_str = f"({H1_r_comp_val}) * î_r + ({H1_theta_comp_val}) * î_θ"
    
    print(f"   H = {sp.pretty(H1_coeff_val, use_unicode=False)} * [ ({sp.pretty(H1_r_comp_val, use_unicode=False)}) î_r + ({sp.pretty(H1_theta_comp_val, use_unicode=False)}) î_θ ]")
    print("\nWhich is equivalent to the expression in Option B:")
    print(f"    H = M_0 * ((2*R_p**3 + R**3) / (3*R**3)) * (-cos(theta) î_r + sin(theta) î_θ)")

    print("-" * 50)
    
    # --- Region 2: Between shield and conductor (R_p < r < R) ---
    print("In the region R_p < r < R:")

    # H-field vector components in the space between
    H2_r_coeff_term1 = (Rp / R)**3
    H2_r_coeff_term2 = (Rp / r)**3
    H2_r_comp_coeff_val = - (2 * M0 / 3) * (H2_r_coeff_term1 - H2_r_coeff_term2)

    H2_theta_coeff_term1 = 2 * (Rp / R)**3
    H2_theta_coeff_term2 = (Rp / r)**3
    H2_theta_comp_coeff_val = (M0 / 3) * (H2_theta_coeff_term1 + H2_theta_coeff_term2)

    print("The components of the H-field are:")
    # Print each number/term in the final equation
    print(f"Hr = - (2 * M_0 / 3) * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta)")
    print(f"   = - (2/3) * M_0 * [ {sp.pretty(H2_r_coeff_term1, use_unicode=False)} - {sp.pretty(H2_r_coeff_term2, use_unicode=False)} ] * cos(theta)\n")
    
    print(f"Hθ = (M_0 / 3) * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(theta)")
    print(f"   = (1/3) * M_0 * [ {sp.pretty(H2_theta_coeff_term1, use_unicode=False)} + {sp.pretty(H2_theta_coeff_term2, use_unicode=False)} ] * sin(theta)")


display_magnetic_field_solution()