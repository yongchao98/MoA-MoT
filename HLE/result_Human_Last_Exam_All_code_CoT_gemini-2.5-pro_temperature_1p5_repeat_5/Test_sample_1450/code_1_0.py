import scipy.special as sp
import scipy.optimize as opt
import numpy as np

def solve_fugacity():
    """
    Solves for the fugacity z where the number density of a Fermi gas
    is 75% of a classical ideal gas under the same pressure and temperature.
    """

    # The Fermi-Dirac integral f_n(z) is related to the polylogarithm function
    # by f_n(z) = -Li_n(-z). We use scipy.special.polylog for this.
    def fermi_dirac_integral(n, z):
        return -sp.polylog(n, -z)

    # We need to solve the equation: f_{3/2}(z) = 0.75 * f_{5/2}(z)
    # This is equivalent to finding the root of: f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0
    def equation_to_solve(z):
        return fermi_dirac_integral(1.5, z) - 0.75 * fermi_dirac_integral(2.5, z)

    try:
        # We need to find a non-trivial root (z > 0).
        # By testing values, we find the root lies between 2 and 3.
        # equation_to_solve(2) is positive, equation_to_solve(3) is negative.
        solution = opt.root_scalar(equation_to_solve, bracket=[2, 3])
        
        if solution.converged:
            z_val = solution.root
            print(f"The equation to be solved is f_3/2(z) = 0.75 * f_5/2(z).")
            print(f"The value of fugacity 'z' that satisfies this is: {z_val:.3f}\n")

            # Final verification step as requested
            print("To verify, we insert the solution back into the equation.")
            lhs_name = "f_3/2(z)"
            rhs_name = "0.75 * f_5/2(z)"
            
            lhs_val = fermi_dirac_integral(1.5, z_val)
            rhs_val = 0.75 * fermi_dirac_integral(2.5, z_val)

            print(f"For the solved z = {z_val:.3f}, the final equation contains these numbers:")
            print(f"  - Left-hand side: {lhs_name} = {lhs_val:.4f}")
            print(f"  - Right-hand side: {rhs_name} = {rhs_val:.4f}")

        else:
            print("Numerical solver did not converge.")

    except Exception as e:
        print(f"An error occurred: {e}")

solve_fugacity()
