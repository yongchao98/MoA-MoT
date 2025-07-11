import scipy.optimize
import scipy.special

def solve_fugacity():
    """
    Calculates the fugacity z for which the number density of a Fermi gas
    is 75% that of a classical ideal gas at the same P and T.
    """

    # The equation to solve is f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0
    # where f_n(z) = -polylog(n, -z)
    def equation_to_solve(z):
        f_3_2 = -scipy.special.polylog(1.5, -z)
        f_5_2 = -scipy.special.polylog(2.5, -z)
        return f_3_2 - 0.75 * f_5_2

    # Find the root of the equation. Based on preliminary checks,
    # the root is expected to be between 10 and 11.
    try:
        solution = scipy.optimize.root_scalar(
            equation_to_solve,
            bracket=[10, 11],
            method='brentq'
        )
        z = solution.root
    
        # Verify the solution
        f_3_2_val = -scipy.special.polylog(1.5, -z)
        f_5_2_val = -scipy.special.polylog(2.5, -z)
        
        print("Solving the equation: f_{3/2}(z) = 0.75 * f_{5/2}(z)")
        print(f"Found fugacity z = {z:.4f}")
        print(f"At this fugacity:")
        print(f"  f_{{3/2}}({z:.2f}) = {f_3_2_val:.4f}")
        print(f"  0.75 * f_{{5/2}}({z:.2f}) = {0.75 * f_5_2_val:.4f}")

        # Round the result to two significant digits
        # For a number like 10.195, rounding to 2 s.f. gives 10.
        z_rounded = float(f'{z:.1e}')
        if z_rounded == 10.0:
            z_rounded = int(z_rounded) # Display as 10 instead of 10.0

        print(f"\nThe value of the fugacity to two significant digits is: {z_rounded}")

    except (ImportError, ValueError) as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have scipy installed (`pip install scipy`).")
        print("The bracket [10, 11] might not contain the root.")

solve_fugacity()
