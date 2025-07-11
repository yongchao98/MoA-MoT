import mpmath

def solve_fermi_fugacity():
    """
    Calculates the fugacity z for a Fermi gas under the specified conditions.

    The condition is that the number density of the Fermi gas is 75% that of a
    classical ideal gas under the same pressure and temperature. This leads to
    the equation: f_{3/2}(z) = 0.75 * f_{5/2}(z), where f_n(z) is the
    Fermi-Dirac integral of order n.
    """
    # Set precision for mpmath calculations for accuracy
    mpmath.mp.dps = 25

    # The Fermi-Dirac integral f_n(z) is related to the polylogarithm function Li_n(x)
    # by f_n(z) = -Li_n(-z). In mpmath, Li_n(x) is mpmath.polylog(n, x).
    def fermi_dirac_integral(n, z):
        """Computes the Fermi-Dirac integral f_n(z)."""
        return -mpmath.polylog(n, -z)

    # We need to solve f_{3/2}(z) = 0.75 * f_{5/2}(z).
    # This is equivalent to finding the root of g(z) = f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0.
    
    # Define the numbers in the equation
    n1 = 1.5
    n2 = 2.5
    C = 0.75

    def equation_to_solve(z):
        f_n1 = fermi_dirac_integral(n1, z)
        f_n2 = fermi_dirac_integral(n2, z)
        return f_n1 - C * f_n2

    print(f"The equation to be solved for the fugacity (z) is:")
    print(f"f_{n1}(z) = {C} * f_{n2}(z)\n")

    # Use a numerical root finder. An initial guess of z=1.0 is reasonable.
    fugacity_solution = mpmath.findroot(equation_to_solve, 1.0)
    
    # Round the result to two significant digits.
    num_significant_digits = 2
    rounded_fugacity = mpmath.nstr(fugacity_solution, num_significant_digits)

    print(f"The numerically calculated value of the fugacity is: {fugacity_solution}")
    print(f"The fugacity rounded to two significant digits is: {rounded_fugacity}")

solve_fermi_fugacity()