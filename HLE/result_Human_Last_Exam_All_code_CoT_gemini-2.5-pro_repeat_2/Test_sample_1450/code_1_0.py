import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import gamma

def fermi_dirac_integral(n, z):
    """
    Computes the complete Fermi-Dirac integral f_n(z) using numerical integration.
    f_n(z) = (1/Gamma(n)) * integral from 0 to inf of x^(n-1) / (exp(x)/z + 1) dx
    """
    if z <= 0:
        # Fugacity must be positive.
        return 0
    
    # The integrand function
    integrand = lambda x: np.power(x, n - 1) / (np.exp(x) / z + 1)
    
    # Numerically integrate from 0 to infinity
    integral_val, _ = quad(integrand, 0, np.inf)
    
    return integral_val / gamma(n)

def target_equation(z):
    """
    Defines the equation to be solved: f_{3/2}(z) / f_{5/2}(z) - 0.75 = 0.
    """
    f_3_2 = fermi_dirac_integral(1.5, z)
    f_5_2 = fermi_dirac_integral(2.5, z)
    
    if f_5_2 == 0:
        # This case is for z=0. The ratio is 1, so the function value is 1 - 0.75 = 0.25
        return 0.25
        
    return f_3_2 / f_5_2 - 0.75

try:
    # We solve for the root in the interval [0.1, 10.0]. The ratio of the integrals
    # is a monotonically decreasing function of z from 1 to 0, so a unique
    # solution is guaranteed.
    z_solution = brentq(target_equation, 0.1, 10.0)

    # Now, we display the numbers that make up the final equation.
    f32_val = fermi_dirac_integral(1.5, z_solution)
    f52_val = fermi_dirac_integral(2.5, z_solution)
    ratio = f32_val / f52_val

    print("The problem reduces to solving the equation: f_{3/2}(z) / f_{5/2}(z) = 0.75")
    print(f"A numerical solver found the fugacity z to be: {z_solution}")
    print("\nVerifying the solution by showing the values in the equation:")
    print(f"For z = {z_solution:.5f}:")
    print(f"  f_3/2(z) = {f32_val:.5f}")
    print(f"  f_5/2(z) = {f52_val:.5f}")
    print(f"  Ratio = {f32_val:.5f} / {f52_val:.5f} = {ratio:.5f}")
    
    # The final answer must be given to two significant digits.
    z_rounded = float(f"{z_solution:.2g}")
    print(f"\nThe value of the fugacity z to two significant digits is: {z_rounded}")

except Exception as e:
    print(f"An error occurred during the calculation: {e}")

<<<1.4>>>