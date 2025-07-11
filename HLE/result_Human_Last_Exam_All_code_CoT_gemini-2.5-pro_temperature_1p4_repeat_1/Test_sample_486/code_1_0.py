import numpy as np
from scipy.integrate import quad

def calculate_integral_asymptotics(R_values):
    """
    Calculates the value of I(R)/R^2 for the 1D solution for a list of radii R.
    I(R) = integral_{B_R} |grad u|^2 dV
    For the 1D solution u(x1) = tanh(x1/sqrt(2)), this becomes:
    I(R) = pi * integral_{-R to R} (R^2 - x^2) * 0.5 * sech^4(x/sqrt(2)) dx
    """
    
    # The integrand function for a given R.
    # Note: sech(z) = 1/cosh(z)
    def integrand(x, R):
        sech_val = 1.0 / np.cosh(x / np.sqrt(2))
        return np.pi * (R**2 - x**2) * 0.5 * sech_val**4

    print("Demonstrating the convergence of R^{-2} * integral(|\nabla u|^2) for the 1D solution.")
    print("-" * 70)
    print(f"{'Radius R':<10} | {'Integral I(R)':<25} | {'I(R) / R^2':<25}")
    print("-" * 70)

    for R in R_values:
        # Perform numerical integration from -R to R
        integral_value, error = quad(integrand, -R, R, args=(R,))
        
        # Calculate the ratio
        ratio = integral_value / (R**2)
        
        print(f"{R:<10} | {integral_value:<25.10f} | {ratio:<25.10f}")

    # The theoretical limit as R -> infinity
    theoretical_limit = (2 * np.sqrt(2) * np.pi) / 3
    print("-" * 70)
    print(f"Theoretical limit as R -> infinity: {theoretical_limit:.10f}")
    
    # The final answer for a
    a = 2
    print(f"\nThe largest possible value for 'a' is {a}.")

if __name__ == '__main__':
    # Use a set of increasing radii to show the convergence
    radii = [10, 50, 100, 200, 500]
    calculate_integral_asymptotics(radii)
