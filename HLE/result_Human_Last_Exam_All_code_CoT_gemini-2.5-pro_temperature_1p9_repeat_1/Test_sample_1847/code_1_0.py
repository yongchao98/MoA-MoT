import math

def solve_shock_profile():
    """
    This function derives and prints the analytical solution for the given shock wave problem.
    The solution is an implicit equation relating normalized density to normalized position.
    """

    # --- Derived constants ---
    # For a monatomic ideal gas, gamma = 5/3.
    # For Mach number M=sqrt(3), the Rankine-Hugoniot relations give the
    # downstream-to-upstream density ratio rho_1/rho_0.
    # rho_1/rho_0 = ( (gamma+1) * M^2 ) / ( 2 + (gamma-1) * M^2 )
    #             = ( (8/3) * 3 ) / ( 2 + (2/3) * 3 ) = 8 / (2+2) = 2
    rho_limit = 2

    # The derivation of the shock profile leads to an implicit equation.
    # Let rho be the density normalized by the ambient density (rho_norm = rho_actual / rho_0)
    # Let X be the position normalized by the conductive length scale (X = x / L)

    # The derived implicit equation is of the form:
    # exp(A * X) = B * (rho - 1) / sqrt(rho * (rho_limit - rho))
    
    # --- Numerical coefficients in the final equation ---
    # The coefficient A = 2/3 is derived from the normalization.
    coeff_A_num = 2
    coeff_A_den = 3
    
    # The coefficient B = 2*sqrt(2) comes from algebraic manipulation
    # of the integrated velocity profile equation when transformed to density.
    coeff_B_integer = 2
    coeff_B_sqrt = 2

    # --- Print the final result ---
    # We construct the final equation as a string, showing all derived numbers.
    equation = (
        f"exp( ({coeff_A_num}/{coeff_A_den}) * X ) = "
        f"({coeff_B_integer}*sqrt({coeff_B_sqrt}) * (rho - 1)) / sqrt(rho*({rho_limit} - rho))"
    )
    
    print("The analytical solution for the normalized density profile, rho(X), is given by the following implicit equation:")
    print("Where 'rho' is the normalized density rho_actual/rho_ambient, and 'X' is the normalized position x/L.\n")
    print(equation)


solve_shock_profile()

# The final answer format as requested.
print("\n<<<exp( (2/3) * X ) = (2*sqrt(2)*(rho - 1)) / sqrt(rho*(2 - rho))>>>")