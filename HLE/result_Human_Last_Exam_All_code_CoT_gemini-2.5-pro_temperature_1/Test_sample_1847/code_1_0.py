import sympy

def solve_shock_profile():
    """
    This function derives and prints the analytical solution for the density profile
    of a specific shock wave.

    The problem asks for the density profile of a Mach sqrt(3) shock wave in a
    monatomic ideal gas (gamma = 5/3) with a Prandtl number of 3/4.
    
    The analytical solution for this case, known as Becker's solution, is an
    implicit equation relating the non-dimensional position (x') to the
    non-dimensional density (rho').

    - rho' = rho / rho_0 (density normalized by the upstream ambient density)
    - x' = x / L, where L is the given length scale.
      L = kappa / (rho_0 * M * c_0 * C_v)
    
    The equation defines the shape of the shock wave, with rho' transitioning
    from 1 (upstream) to 2 (downstream) as x' goes from -infinity to +infinity.
    """

    # Define symbols for clarity in the output
    x_prime = sympy.Symbol("x/L")
    rho_prime = sympy.Symbol("rho/rho_0")

    # The pre-factor in the solution is calculated as (8 * Pr) / (3 * (gamma + 1))
    # Given Pr = 3/4 and gamma = 5/3
    # pre_factor = (8 * (3/4)) / (3 * (5/3 + 1)) = 6 / (3 * 8/3) = 6 / 8 = 3/4
    pre_factor_num = 3
    pre_factor_den = 4
    pre_factor = sympy.Rational(pre_factor_num, pre_factor_den)

    # The argument of the logarithm is derived from integrating the momentum equation
    # after the simplification from the energy equation.
    # The numbers 8, 2, and 1 arise from the specific shock conditions (M=sqrt(3))
    # and the choice of the shock center definition.
    numerator_val = 8
    num_term_val_1 = 1
    den_term_val_1 = 2
    
    numerator = numerator_val * (rho_prime - num_term_val_1)**2
    denominator = rho_prime * (den_term_val_1 - rho_prime)
    log_argument = numerator / denominator

    # Construct the final equation: x' = pre_factor * ln(log_argument)
    final_equation = sympy.Eq(x_prime, pre_factor * sympy.log(log_argument))

    # Print the explanation and the final equation
    print("The analytical solution for the density profile is an implicit equation relating")
    print("the non-dimensional position (x/L) and the non-dimensional density (rho/rho_0).")
    print("\nFor a Mach sqrt(3) shock in a monatomic gas with Pr = 3/4, the equation is:\n")
    
    # Use pretty print for a more readable equation format
    sympy.pprint(final_equation, use_unicode=False)

    print("\nWhere:")
    print("  - rho is the local density.")
    print("  - rho_0 is the ambient density far upstream of the shock.")
    print("  - x is the position coordinate.")
    print("  - L is the ambient conductive length scale L = kappa/(rho_0 * M * c_0 * C_v).")

if __name__ == '__main__':
    solve_shock_profile()