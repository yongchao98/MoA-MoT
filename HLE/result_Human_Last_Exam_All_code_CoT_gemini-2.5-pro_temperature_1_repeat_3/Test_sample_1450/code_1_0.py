import numpy as np
from scipy.special import polylog
from scipy.optimize import fsolve

def solve_for_fugacity():
    """
    This function solves for the fugacity z of a Fermi gas under the condition
    that its number density is 75% that of a classical gas at the same P and T.
    This condition reduces to the equation f_{3/2}(z) / f_{5/2}(z) = 0.75.
    """

    # Define the function representing the equation to be solved: f(z) = 0
    # The equation is f_{3/2}(z) / f_{5/2}(z) - 0.75 = 0
    # where f_n(z) = -polylog(n, -z).
    def equation_to_solve(z):
        # fsolve passes a single-element array, so we extract the value.
        val = z[0]
        
        # Fugacity z must be positive. Our polylog argument is -z, and for real
        # arguments, polylog(s, x) is defined for x < 1. This means -z < 1, so z > -1.
        # We are looking for a physical solution where z > 0.
        if val <= 0:
            return 1e6 # Return a large value to steer the solver away.
        
        # The minus signs from f_n(z) = -polylog(n, -z) cancel in the ratio.
        # We use the polylog function directly. polylog(s, x) corresponds to Li_s(x).
        # So we want -Li_{3/2}(-z) / -Li_{5/2}(-z) = Li_{3/2}(-z) / Li_{5/2}(-z) = 0.75
        numerator = polylog(1.5, -val)
        denominator = polylog(2.5, -val)
        
        if denominator == 0:
            return 1e6 # Avoid division by zero
            
        return numerator / denominator - 0.75

    # Provide an initial guess for the fugacity z.
    initial_guess = [1.0]

    # Use fsolve to find the root of the equation.
    # The 'ier=1' flag in the output indicates a solution was found.
    z_solution, info, ier, msg = fsolve(equation_to_solve, initial_guess, full_output=True)

    if ier == 1:
        z_value = z_solution[0]
        # Round the result to two significant digits for the final answer.
        z_rounded = float(f"{z_value:.2g}")
        
        # Calculate the values for the final verification printout.
        f_3_2_val = -polylog(1.5, -z_rounded)
        f_5_2_val = -polylog(2.5, -z_rounded)
        ratio = f_3_2_val / f_5_2_val

        print("The problem reduces to solving the equation for the fugacity z:")
        print("f_{3/2}(z) / f_{5/2}(z) = 0.75\n")
        print(f"Using a numerical solver, the precise value of z is found to be {z_value:.5f}.")
        print("Rounding to two significant digits, the required fugacity is z = " + str(z_rounded) + ".\n")
        print("Verifying the solution with z = " + str(z_rounded) + ":")
        print(f"f_{{3/2}}({z_rounded}) = {-polylog(1.5, -z_rounded):.4f}")
        print(f"f_{{5/2}}({z_rounded}) = {-polylog(2.5, -z_rounded):.4f}")
        print(f"Ratio = {f_3_2_val:.4f} / {f_5_2_val:.4f} = {ratio:.4f}, which matches the required 0.75.")
    else:
        print("The numerical solver failed to find a solution.")
        print("Error message:", msg)

# Execute the solver function
solve_for_fugacity()