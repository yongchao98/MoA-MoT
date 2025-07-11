import numpy as np
from scipy.special import polylog
from scipy.optimize import root_scalar

def equation_to_solve(z):
    """
    This function represents the equation we want to solve:
    f_3/2(z) - 0.75 * f_5/2(z) = 0
    where the Fermi-Dirac integral f_n(z) is given by -polylog(n, -z).
    """
    # f_3/2(z) is -polylog(1.5, -z)
    term1 = -polylog(1.5, -z)
    # f_5/2(z) is -polylog(2.5, -z)
    term2 = -polylog(2.5, -z)
    
    # The equation is term1 - 0.75 * term2 = 0
    return term1 - 0.75 * term2

# We need to solve the equation for z.
# Based on asymptotic analysis, we know the root is a positive number larger than 1.
# We can use a numerical root finder with a suitable bracket, e.g., [1, 100].
try:
    solution = root_scalar(equation_to_solve, bracket=[1, 100])
    z = solution.root
    
    # Format the result to two significant digits.
    # The 'g' format specifier is suitable for this.
    z_rounded_str = f"{z:.2g}"

    print("For an ideal Fermi gas, the number density (ν) and pressure (P) are:")
    print("ν ∝ f_3/2(z)")
    print("P/(k_B T) ∝ f_5/2(z)")
    print("\nFor a classical ideal gas, ν_cl = P/(k_B T).")
    print("The condition ν = 0.75 * ν_cl thus implies the following equation for the fugacity z:")
    print("f_3/2(z) = 0.75 * f_5/2(z)")
    print("\nSolving this equation numerically yields:")
    print(f"z = {z:.4f}")
    print(f"\nThe value of the fugacity z to two significant digits is: {z_rounded_str}")

except (ImportError, ModuleNotFoundError):
    print("This script requires the 'scipy' library.")
    print("Please install it using: pip install scipy")
except Exception as e:
    print(f"An error occurred: {e}")
