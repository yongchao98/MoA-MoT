import numpy as np
from scipy.special import polylog
from scipy.optimize import root_scalar

def solve_for_fugacity():
    """
    This function finds the fugacity z for a Fermi gas under the condition that
    its number density is 75% that of a classical ideal gas at the same P and T.
    
    This reduces to solving the equation: f_{3/2}(z) / f_{5/2}(z) = 0.75
    where f_n(z) is the Fermi-Dirac integral, which is related to the polylogarithm
    function Li_n by f_n(z) = -polylog(n, -z).
    """

    # Define the function whose root we want to find.
    # We are solving: f_{3/2}(z) / f_{5/2}(z) - 0.75 = 0
    # which is equivalent to: polylog(1.5, -z) / polylog(2.5, -z) - 0.75 = 0
    def equation_to_solve(z):
        # The fugacity z must be positive.
        if z <= 0:
            return np.inf # Return a large number for invalid domain
        
        # Use the polylog function from scipy.special.
        # Note: f_n(z) = -polylog(n, -z)
        # The ratio of two f_n(z) will have the negative signs cancel out.
        f_3_2 = -polylog(1.5, -z)
        f_5_2 = -polylog(2.5, -z)
        
        # Avoid division by zero, although f_{5/2}(z) is non-zero for z > 0
        if f_5_2 == 0:
            return np.inf
            
        return (f_3_2 / f_5_2) - 0.75

    # We need to find the root of the function. We can determine a bracket for the root.
    # At z -> 0, the ratio is 1, so equation_to_solve(z) -> 0.25 > 0.
    # For large z, the ratio is approx. 2.5/ln(z).
    # 2.5/ln(z) = 0.75 => ln(z) = 10/3 => z is approx. e^(10/3) ~= 28.
    # Let's check a value larger than this, e.g., z=30.
    # equation_to_solve(30) is negative.
    # So, a root exists in the bracket [1, 30].
    try:
        sol = root_scalar(equation_to_solve, bracket=[1, 30], method='brentq')
        fugacity_z = sol.root

        # Print the final equation as requested
        print("The problem is solved by finding the fugacity z in the following equation:")
        print("f_{3/2}(z) / f_{5/2}(z) = 0.75")
        print("\nwhere f_n(z) is the Fermi-Dirac integral.")
        
        # Print the result formatted to two significant digits
        print(f"\nThe value of the fugacity z for which this is true is: {fugacity_z:.2g}")

    except Exception as e:
        print(f"An error occurred during the numerical solution: {e}")

# Run the solver
solve_for_fugacity()