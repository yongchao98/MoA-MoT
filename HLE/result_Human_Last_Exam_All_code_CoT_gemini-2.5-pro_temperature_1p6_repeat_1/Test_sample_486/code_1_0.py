import numpy as np
from scipy.integrate import quad

def analyze_integral_growth():
    """
    Numerically verifies the asymptotic growth of the Dirichlet integral
    for the solution u(x,y,z) = tanh(z/sqrt(2)).
    """
    
    # The integrand is pi*(R^2 - z^2) * (1/2) * sech^4(z/sqrt(2))
    # We define the part without the constant factors for the numerical integration.
    def integrand(z, R):
        # sech(x) = 1 / cosh(x)
        sech_val = 1 / np.cosh(z / np.sqrt(2))
        return (R**2 - z**2) * (sech_val**4)

    def I(R):
        # Full integral value. Factor is (pi/2).
        integral_val, _ = quad(integrand, -R, R, args=(R,))
        return (np.pi / 2) * integral_val

    # Use a range of R values that are large enough to see the asymptotic behavior.
    R_values = np.logspace(1, 2, 8)  # R from 10 to 100

    # Calculate the integral for each R.
    I_values = np.array([I(R) for R in R_values])

    # To find the exponent 'a' in I(R) ≈ C*R^a, we can fit a line to
    # log(I(R)) vs log(R). The slope of the line will be 'a'.
    log_R = np.log(R_values)
    log_I = np.log(I_values)
    
    # Perform a linear fit (degree 1 polynomial)
    coeffs = np.polyfit(log_R, log_I, 1)
    a_numerical = coeffs[0]

    print("This script numerically verifies the scaling of the integral I(R).")
    print("The integral I(R) is expected to scale as C * R^a.")
    print(f"By fitting a line to log(I) vs log(R), we estimate 'a'.\n")
    print(f"Numerical estimate for the exponent 'a': {a_numerical:.6f}\n")

    # The asymptotic equation from the analytical derivation is:
    # Integral(R) ≈ (2 * pi * sqrt(2) / 3) * R^2
    print("The analytical derivation gives the asymptotic equation:")
    print("Integral(R) \u2248 C * R^a, where C = (2 * \u03c0 * \u221a2) / 3 and a = 2.\n")
    
    print("Printing the numbers from the final analytical equation:")
    print(f"Coefficient numerator factor 1: 2")
    print(f"Coefficient numerator factor 2 (pi): {np.pi}")
    print(f"Coefficient numerator factor 3 (sqrt(2)): {np.sqrt(2)}")
    print(f"Coefficient denominator: 3")
    print(f"Exponent 'a': 2")

analyze_integral_growth()