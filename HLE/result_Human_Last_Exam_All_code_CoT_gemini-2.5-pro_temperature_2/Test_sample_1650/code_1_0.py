import math

def display_overlap_integral_2s():
    """
    This function displays the analytical expression for the overlap integral
    of two 2s orbitals in a diatomic molecule like H2+.
    """
    
    # Define the components of the overlap integral expression
    # based on the analytical derivation.
    # The result is in terms of the internuclear distance R and
    # the effective nuclear charge zeta.
    
    # Let rho = zeta * R / 2
    # The expression for S in terms of rho is:
    # S = exp(-rho) * (1 + rho + rho^2/3 + rho^4/15)
    
    # Substituting rho = zeta*R/2 back into the equation:
    term0 = "exp(-zeta*R/2)"
    term1_coeff = 1
    term2_coeff_num = 1
    term2_coeff_den = 2
    term3_coeff_num = 1
    term3_coeff_den = 12
    term4_coeff_num = 1
    term4_coeff_den = 240

    print("The overlap integral S for two 2s orbitals is:")
    print("S(R, zeta) = ", end="")
    print(f"{term0} * (", end="")
    print(f"{term1_coeff} ", end="")
    print(f"+ ({term2_coeff_num}/{term2_coeff_den})*zeta*R ", end="")
    print(f"+ ({term3_coeff_num}/{term3_coeff_den})*zeta^2*R^2 ", end="")
    print(f"+ ({term4_coeff_num}/{term4_coeff_den})*zeta^4*R^4", end="")
    print(")")
    print("\nWhich simplifies to:")
    final_expression = (f"S(R, zeta) = exp(-zeta*R/2) * (1 + (zeta*R)/2 + "
                        f"(zeta^2*R^2)/{term3_coeff_den} + (zeta^4*R^4)/{term4_coeff_den})")
    print(final_expression)

display_overlap_integral_2s()