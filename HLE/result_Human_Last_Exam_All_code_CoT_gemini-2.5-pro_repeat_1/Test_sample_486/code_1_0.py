import numpy as np
from scipy.integrate import quad

def solve_and_verify():
    """
    Numerically verifies the growth rate of the Dirichlet energy for a 1D solution
    to the given PDE, and determines the value of 'a'.
    """
    print("This script analyzes the growth rate of the integral of |nabla u|^2.")
    print("We consider the 1D solution u(x,y,z) = tanh(z/sqrt(2)).")
    print("For this solution, |nabla u|^2 = 0.5 * sech(z/sqrt(2))^4.")
    print("The integral over a ball B_R is given by:")
    print("I(R) = Integral from -R to R of [pi * (R^2 - z^2) * 0.5 * sech(z/sqrt(2))^4] dz\n")

    # The integrand function for the 1D integral in z
    # I(R) = integral_{-R to R} integrand(z, R) dz
    def integrand(z, R):
        return np.pi * (R**2 - z**2) * 0.5 * (1 / np.cosh(z / np.sqrt(2)))**4

    radii = [10.0, 20.0, 50.0, 100.0, 200.0]

    print("Numerically calculating the integral I(R) for increasing R:")
    print("-" * 50)
    print(f"{'R':>10s} {'Integral I(R)':>20s} {'I(R) / R^2':>15s}")
    print("-" * 50)

    for R in radii:
        # Perform the numerical integration using scipy.integrate.quad
        integral_value, error = quad(integrand, -R, R, args=(R,))
        ratio = integral_value / (R**2)
        print(f"{R:>10.1f} {integral_value:>20.4f} {ratio:>15.4f}")

    # Theoretical limit of the ratio I(R)/R^2 as R -> infinity
    # The limit is (pi/2) * integral_{-inf to inf} sech^4(z/sqrt(2)) dz
    # integral_{-inf to inf} sech^4(y) dy = 4/3
    # integral_{-inf to inf} sech^4(z/sqrt(2)) dz = sqrt(2) * 4/3
    # Limit = (pi/2) * sqrt(2) * 4/3 = 2*pi*sqrt(2)/3
    theoretical_limit = 2 * np.pi * np.sqrt(2) / 3
    print("-" * 50)
    print(f"The theoretical limit of I(R)/R^2 is: {theoretical_limit:.4f}\n")
    
    print("The numerical results show that the integral grows proportionally to R^2.")
    print("This means the ratio R^(-a) * I(R) converges to a positive, finite number if a = 2.")
    print("For the limit inferior to be > 0, we need a <= 2.")
    
    a = 2
    
    # Per instructions, outputting the final equation and answer
    print("\nThe problem is to find the largest 'a' such that:")
    print(f"liminf_{{R->inf}} R**(-a) * Integral_B_R(|nabla u|^2) > 0")
    print(f"\nThe largest possible value for a is: {a}")


solve_and_verify()
<<<2>>>