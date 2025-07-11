import numpy as np
from scipy.special import polylog
from scipy.optimize import root_scalar

def solve_fugacity():
    """
    This function sets up and solves the equation for the fugacity z and prints the result.
    """
    print("Step 1: The relationship between the number density of a Fermi gas (ν) and a classical ideal gas (ν_cl) is:")
    print("ν / ν_cl = f_3/2(z) / f_5/2(z)")
    print("where z is the fugacity and f_n(z) is the Fermi-Dirac integral.")
    print("-" * 40)

    print("Step 2: From the problem statement (ν = 75% of ν_cl), we derive the equation to solve for z:")
    # The numbers in the final equation are 1.5 (or 3/2), 0.75, and 2.5 (or 5/2).
    n1 = 1.5
    n2 = 2.5
    factor = 0.75
    print(f"f_{n1}(z) - {factor} * f_{n2}(z) = 0")
    print("-" * 40)

    def equation_to_solve(z: float, n1: float, n2: float, factor: float) -> float:
        """
        Defines the function f_n1(z) - factor * f_n2(z) whose root is the solution.
        It uses the relation f_n(z) = -polylog(n, -z).
        """
        # Fugacity z must be a positive real number.
        if z <= 0:
            return np.inf

        # Calculate the two terms of the equation using the polylogarithm function
        term1 = -polylog(n1, -z)
        term2 = -factor * (-polylog(n2, -z))
        
        # We seek a real root, so we return the real part of the sum.
        return np.real(term1 + term2)

    print("Step 3: Solving the equation numerically using SciPy's root-finding capabilities.")
    # A bracket is chosen based on preliminary analysis showing the root is between 1 and 2.
    try:
        solution = root_scalar(equation_to_solve, args=(n1, n2, factor), bracket=[0.1, 5.0], method='brentq')
        result_z = solution.root
        print(f"The numerical solver found a solution at z = {result_z:.5f}")
    except (ValueError, RuntimeError) as e:
        print(f"The solver failed to find a solution: {e}")
        return

    print("-" * 40)
    print("Step 4: Reporting the final answer.")
    # Format the result to two significant digits as requested.
    # The '.2g' format specifier handles this correctly (e.g., 1.86... -> 1.9).
    final_answer = float(f"{result_z:.2g}")
    print(f"The value of the fugacity, z, to two significant digits is: {final_answer}")
    print(f"<<<{final_answer}>>>")

if __name__ == "__main__":
    solve_fugacity()
