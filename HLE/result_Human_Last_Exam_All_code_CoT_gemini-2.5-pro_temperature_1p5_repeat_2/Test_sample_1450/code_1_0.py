import numpy as np
from scipy.integrate import quad
from scipy.optimize import root_scalar
from scipy.special import gamma

def fermi_dirac_integral(n, z):
    """
    Computes the Fermi-Dirac integral f_n(z) using numerical integration.
    f_n(z) = (1/Gamma(n)) * integral from 0 to inf of (x^(n-1)) / (exp(x)/z + 1) dx.
    This is also written as -Li_n(-z).
    """
    if z <= 0:
        return 0.0
    
    # The integrand is x^(n-1) / (exp(x - log(z)) + 1).
    # Using log(z) improves numerical stability for a wide range of z.
    log_z = np.log(z)
    integrand = lambda x: (x**(n - 1)) / (np.exp(x - log_z) + 1)
    
    # Perform numerical integration from 0 to infinity.
    result, _ = quad(integrand, 0, np.inf)
    
    # The integral is scaled by 1/Gamma(n).
    return result / gamma(n)

def target_function_to_solve(z):
    """
    Returns the value of f_{3/2}(z) - 0.75 * f_{5/2}(z).
    The root of this function is the solution for the fugacity z.
    """
    f_3_2 = fermi_dirac_integral(1.5, z)
    f_5_2 = fermi_dirac_integral(2.5, z)
    return f_3_2 - 0.75 * f_5_2

# Find the root of the target function. A root must exist between z where the
# function is positive (e.g., small z) and where it is negative (large z).
# A bracket of [0.1, 5] is sufficient.
try:
    solution = root_scalar(target_function_to_solve, bracket=[0.1, 5], method='brentq')
    z_value = solution.root

    # Calculate the final values for printing verification
    f_3_2_final = fermi_dirac_integral(1.5, z_value)
    f_5_2_final = fermi_dirac_integral(2.5, z_value)

    print("The problem requires finding the fugacity 'z' such that the number density of a Fermi gas")
    print("is 75% that of a classical gas under the same pressure P and temperature T.\n")
    print("This condition leads to the equation relating Fermi-Dirac integrals:")
    print("ν_fermi / ν_classical = f_{3/2}(z) / f_{5/2}(z) = 0.75\n")

    print(f"Solving for z numerically gives a fugacity value of: {z_value:.2g}\n")

    print("Verification using the solved value of z:")
    print(f"For z = {z_value:.4f}:")
    print(f"  f_3/2({z_value:.4f}) = {f_3_2_final:.4f}")
    print(f"  f_5/2({z_value:.4f}) = {f_5_2_final:.4f}")
    print(f"Their ratio is {f_3_2_final:.4f} / {f_5_2_final:.4f} = {f_3_2_final / f_5_2_final:.4f}, which is approximately 0.75.")

except (ImportError, ModuleNotFoundError):
    print("This script requires the 'scipy' and 'numpy' libraries.")
    print("Please install them using: pip install scipy numpy")
except ValueError as e:
    print(f"An error occurred during calculation: {e}")
    print("This may happen if the root is not within the specified bracket.")