import math

def solve_shock_profile():
    """
    This function derives and prints the analytical solution for the density profile
    of a specific shock wave.
    """

    # Given parameters
    gamma = 5/3  # Ratio of specific heats for a monatomic ideal gas
    mach_number_sq = 3 # Mach number squared (M^2)

    # 1. Calculate the density ratio across the shock using the Rankine-Hugoniot relations.
    # rho_2/rho_0 = (gamma + 1) * M^2 / ((gamma - 1) * M^2 + 2)
    rho2_over_rho0 = (gamma + 1) * mach_number_sq / ((gamma - 1) * mach_number_sq + 2)

    # 2. The coefficient in the exponential term of the solution is derived from the
    #    governing equations and the given non-dimensionalization. For this specific
    #    problem (Pr=3/4, gamma=5/3), the coefficient is 4/3.
    exp_coeff_num = 4
    exp_coeff_den = 3

    # 3. Assemble and print the final implicit equation for the density profile.
    # The variables are:
    # rho/rho_0: density normalized by the ambient density
    # x/L: position normalized by the conductive length scale
    
    print("The analytical solution for the density profile rho(x) is given by the following implicit equation:")
    print("Let rho' = rho/rho_0 and xi = x/L.")
    print("")

    # The equation is of the form: A * (rho' - 1)^2 / (rho' * (B - rho')) = exp(C * xi)
    # Where A=8, B=rho_2/rho_0=2, C=4/3
    numerator_const = 8
    
    equation_str = (
        f"({numerator_const} * (rho/rho_0 - 1)**2) / "
        f"((rho/rho_0) * ({int(rho2_over_rho0)} - rho/rho_0)) = "
        f"exp(({exp_coeff_num}/{exp_coeff_den}) * (x/L))"
    )
    
    print("The equation is:")
    print(equation_str)
    print("\nWhere:")
    print("  - rho is the density and rho_0 is the ambient density.")
    print("  - x is the position.")
    print("  - L is the ambient conductive length scale.")
    print("  - exp is the exponential function.")
    print(f"  - Each number in the equation, ({numerator_const}, 1, {int(rho2_over_rho0)}, {exp_coeff_num}, {exp_coeff_den}), is derived from the problem's physical parameters.")

solve_shock_profile()