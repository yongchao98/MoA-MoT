import scipy.optimize
import math

def solve_fugacity():
    """
    Solves for the fugacity z where the number density of a Fermi gas is
    75% of a classical ideal gas under the same pressure and temperature.
    """
    
    # The problem reduces to solving the equation:
    # f_3/2(z) - 0.75 * f_5/2(z) = 0
    #
    # Using the series expansion for the Fermi-Dirac integrals f_n(z),
    # f_n(z) = sum_{l=1 to inf} [(-1)^(l-1) * z^l / l^n],
    # the equation becomes:
    # sum_{l=1 to inf} [(-1)^(l-1) * z^l * (1/l^(3/2) - 0.75/l^(5/2))] = 0
    #
    # For a non-trivial solution (z != 0), we can divide by z, resulting in an
    # equation F(z) = 0, where F(z) is defined below.

    def equation_to_solve(z, num_terms=100):
        """
        Defines the function F(z) whose root we need to find.
        F(z) = sum_{l=1 to inf} [(-1)^(l-1) * z^(l-1) * ((l - 0.75) / l^(5/2))]
        """
        total = 0.0
        for l in range(1, num_terms + 1):
            term = ((-1)**(l - 1)) * (z**(l - 1)) * ((l - 0.75) / (l**2.5))
            total += term
        return total

    # Find the root of the function. Based on preliminary analysis, the root
    # is expected to be between 1.0 and 1.5.
    try:
        sol = scipy.optimize.root_scalar(equation_to_solve, bracket=[1.0, 1.5], method='brentq')
        z_solution = sol.root

        # Print the final results as requested
        print("The problem requires solving the equation for the fugacity z:")
        print("ν_fermi = 0.75 * ν_classical")
        print("This simplifies to:")
        print("f_3/2(z) = 0.75 * f_5/2(z)")
        print("\nSolving this equation numerically gives the value of fugacity.")
        # We want to show the final equation with the numbers
        # To format the answer to two significant digits, we can use "{:.2g}".
        print(f"The resulting value for the fugacity z is: {z_solution:.2g}")
        
    except (ImportError, ValueError) as e:
        print("Error: Could not solve the equation.", e)
        print("Please ensure the 'scipy' library is installed ('pip install scipy').")


solve_fugacity()
