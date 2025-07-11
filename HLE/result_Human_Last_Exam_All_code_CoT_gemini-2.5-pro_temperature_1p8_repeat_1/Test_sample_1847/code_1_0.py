import math

def print_solution():
    """
    This function prints the analytical solution for the specified shock wave density profile.
    """
    
    # The analytical solution for the density profile rho(x) is derived from the
    # Becker solution of the Navier-Stokes equations for the special case of a
    # monatomic gas (gamma = 5/3) with a Prandtl number of 3/4. For a Mach sqrt(3)
    # shock, the solution takes a particularly simple form.

    # The normalized density is rho_star = rho/rho_0
    # The normalized position is x_star = x/L

    # The velocity profile u/u_0 = v(x_star) can be expressed as:
    # v(x_star) = 3/4 - (1/4) * tanh((5/9) * x_star)

    # Since rho * u is constant, rho/rho_0 = u_0/u = 1/v(x_star).
    # Therefore, the density profile is:
    # rho_star = 1 / (3/4 - (1/4) * tanh((5/9) * x_star))
    # rho_star = 4 / (3 - tanh((5/9) * x_star))

    numerator_overall = 4
    denominator_const = 3
    numerator_arg = 5
    denominator_arg = 9
    
    equation_string = (
        f"The analytical solution for the density profile is:\n\n"
        f"       {numerator_overall}\n"
        f"ρ/ρ₀ = --------------------------------------------\n"
        f"       {denominator_const} - tanh( ({numerator_arg}/{denominator_arg}) * (x/L) )"
    )
    
    print(equation_string)

print_solution()

# Final Answer Check
# We are expressing the density `ρ` in units of the ambient density `ρ₀`
# and the position `x` in units of the ambient conductive length scale `L`.
# The final equation is rho/rho_0 = 4 / (3 - tanh((5/9) * x/L)).
final_equation_form = "rho/rho_0 = 4/(3-tanh((5/9)*(x/L)))"