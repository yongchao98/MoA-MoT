import math

def display_analytical_solution():
    """
    This function prints the analytical solution for the specified shock wave density profile.

    The problem asks for the analytical solution for the density profile (rho/rho_0) of a
    Mach sqrt(3) shock wave in a monatomic ideal gas with Pr = 3/4. The position is
    expressed in units of the ambient conductive length scale L.

    The solution is not an explicit function rho(x), but an implicit one that relates
    the normalized density and the normalized position.
    """

    # The normalized density rho_ratio = rho / rho_0, where rho_0 is the ambient density.
    # The normalized position xi = x / L.
    # The density increases from rho_0 upstream to 2*rho_0 downstream.
    # The solution contains a constant of integration, D, which depends on the choice of the origin for xi.

    # Coefficients in the final equation
    power_numerator = 1.5
    power_denominator1 = 0.75
    power_denominator2 = 0.75
    density_limit = 2.0

    equation_string = f"""
The analytical solution for the density profile is given by the following implicit equation:

ln[ (rho_ratio - 1)^{power_numerator} / (rho_ratio^{power_denominator1} * (density_limit - rho_ratio)^{power_denominator2}) ] = xi + D

Where:
- rho_ratio is the density divided by the ambient density (rho/rho_0).
- xi is the position divided by the characteristic length scale (x/L).
- D is a constant of integration that sets the position of the shock wave.
- ln is the natural logarithm.

Let's print the equation with the numerical values:
ln[ (rho_ratio - 1)^({power_numerator}) / (rho_ratio^({power_denominator1}) * ({density_limit} - rho_ratio)^({power_denominator2})) ] = xi + D
"""
    print(equation_string)

if __name__ == "__main__":
    display_analytical_solution()
