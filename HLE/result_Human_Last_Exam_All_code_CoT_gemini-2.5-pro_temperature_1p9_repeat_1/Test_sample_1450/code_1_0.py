import scipy.special as sp
import scipy.optimize as opt

def solve_for_fugacity():
    """
    This function finds the fugacity 'z' for a Fermi gas under the specified conditions.
    """
    
    # We are solving the equation stemming from the condition that the number density
    # of the Fermi gas is 75% that of a classical ideal gas at the same P and T.
    # The condition simplifies to: f_{3/2}(z) = 0.75 * f_{5/2}(z)
    #
    # In terms of the polylogarithm function Li_j(x), this is:
    # Li_{1.5}(-z) - 0.75 * Li_{2.5}(-z) = 0

    j1 = 1.5
    j2 = 2.5
    factor = 0.75
    
    print("The task is to solve the following equation for the fugacity 'z':")
    print(f"polylog({j1}, -z) - {factor} * polylog({j2}, -z) = 0")
    print("-" * 60)

    def equation_to_solve(z):
        """
        Defines the transcendental equation g(z) = 0 for fugacity z.
        scipy.special.polylog(j, x) calculates the polylogarithm function Li_j(x).
        """
        return sp.polylog(j1, -z) - factor * sp.polylog(j2, -z)

    try:
        # For the root-finding algorithm, we need an interval [a, b] where the function
        # changes sign. A brief analysis shows the root lies between 10 and 20.
        search_interval = [10, 20]
        
        # Brent's method is a robust and fast root-finding algorithm.
        solution_z = opt.brentq(equation_to_solve, search_interval[0], search_interval[1])
        
        print(f"Numerical solution found for z: {solution_z}")
        
        # The 'g' format specifier with a precision of 2 is used to format the
        # result to two significant digits.
        result_2_sig_digits = f"{solution_z:.2g}"
        
        print(f"\nThe value of fugacity 'z' rounded to two significant digits is: {result_2_sig_digits}")

    except (ValueError, RuntimeError) as e:
        print(f"Could not find a solution in the interval {search_interval}. Error: {e}")

if __name__ == "__main__":
    solve_for_fugacity()
