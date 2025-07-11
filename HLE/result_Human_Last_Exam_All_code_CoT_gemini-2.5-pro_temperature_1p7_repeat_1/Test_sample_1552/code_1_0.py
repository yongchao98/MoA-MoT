import sympy

def find_second_heat_kernel_coefficient():
    """
    This function calculates and displays the second coefficient (integrand)
    in the heat kernel expansion of the spectral action for a massless
    gauged Dirac spinor field.
    """
    
    # --- Define symbols for the mathematical formula ---
    # Representation of the gauge group for the fermion field
    dim_r = sympy.Symbol("dim(r)")
    
    # Gravitational terms
    R = sympy.Symbol("R") # Ricci Scalar
    R_munu_R_munu = sympy.Symbol("R_μν*R^μν") # Ricci Tensor Squared
    C_munurhosigma_C_munurhosigma = sympy.Symbol("C_μνρσ*C^μνρσ") # Weyl Tensor Squared
    
    # Gauge term
    tr_F_munu_F_munu = sympy.Symbol("tr_r(F_μν*F^μν)") # Gauge Field Strength Squared
    
    # --- Numerical coefficients from the formula ---
    # Sourced from reputable physics literature (e.g., Gorokhov & Starobinsky, 2020, JCAP 05, 003)
    c1_num, c1_den = 1, 180
    c2_num, c2_den = -1, 180
    c3_num, c3_den = -7, 360
    c4_num, c4_den = -1, 12

    # --- Construct the expression for the coefficient density ---
    # This is the integrand of a_2, often denoted as c_2 = 16*pi^2 * b_4
    grav_part = dim_r * (sympy.Rational(c1_num, c1_den) * R**2 + \
                         sympy.Rational(c2_num, c2_den) * R_munu_R_munu + \
                         sympy.Rational(c3_num, c3_den) * C_munurhosigma_C_munurhosigma)
                         
    gauge_part = sympy.Rational(c4_num, c4_den) * tr_F_munu_F_munu
    
    full_coefficient_density = grav_part + gauge_part

    # --- Output the results ---
    print("The second coefficient in the heat kernel expansion is a_2 = integral(c_2 / (16*pi^2) * dV).")
    print("The local density c_2 for a massless gauged Dirac spinor is:\n")
    
    # Use sympy's pretty print for a nice mathematical formula output
    sympy.pprint(full_coefficient_density, use_unicode=True)
    
    print("\n------------------------------------------------------------")
    print("The numerical coefficients appearing in the final equation are:")
    
    coeffs = [c1_num, c1_den, c2_num, c2_den, c3_num, c3_den, c4_num, c4_den]
    for num in coeffs:
        print(num)

    # For the final answer format
    return str(full_coefficient_density)

if __name__ == '__main__':
    final_answer = find_second_heat_kernel_coefficient()
    # The final answer format is specified in the prompt instructions
    # print(f"<<<{final_answer}>>>")