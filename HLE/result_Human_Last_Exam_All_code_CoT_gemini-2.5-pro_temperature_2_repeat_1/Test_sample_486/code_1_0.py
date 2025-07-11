import numpy as np
from scipy.integrate import quad

def solve():
    """
    This script calculates the limiting constant for a 1D solution to the Allen-Cahn
    equation, which helps determine the largest possible value for 'a'.

    The problem is to find the largest 'a' such that:
    liminf_{R->inf} R^{-a} * Integral_{B_R} |nabla u|^2 > 0

    1. Upper Bound: General theory on Allen-Cahn solutions states that the total energy
       scales like R^(n-1). For n=3, this is R^2. This implies integral |nabla u|^2
       scales at most like R^2, so a <= 2.

    2. Attainability: We check if a=2 is possible using the 1D solution
       u(x,y,z) = tanh(z/sqrt(2)).
       For this u, |nabla u|^2 = (1/2) * sech^4(z/sqrt(2)).

    3. The integral becomes: Integral_{B_R} |nabla u|^2 dx. In cylindrical coordinates,
       this integral is (pi/2) * Integral_{-R}^{R} (R^2 - z^2) * sech^4(z/sqrt(2)) dz.

    4. The limit C = lim_{R->inf} R^{-2} * Integral_{B_R} |nabla u|^2 dx is then:
       C = (pi/2) * Integral_{-inf}^{inf} sech^4(z/sqrt(2)) dz.
       If C is a positive constant, then a=2 is attainable.

    We calculate C both analytically and numerically.
    """

    # Part 1: Calculate the constant from its analytical formula.
    # The analytical result of the integral is (pi/2) * (4*sqrt(2)/3).
    # This simplifies to C = (2 * pi * sqrt(2)) / 3.
    # The following code prints out the components of this "final equation" for the limit C.

    p1 = 2.0
    p2 = np.pi
    p3 = np.sqrt(2.0)
    denominator = 3.0

    print("--- Analytical Calculation of the Limit Constant C ---")
    print(f"The analytical expression for the limit is (p1 * pi * sqrt(2)) / denominator.")
    print(f"Component p1: {p1}")
    print(f"Component pi (p2): {p2}")
    print(f"Component sqrt(2) (p3): {p3}")
    print(f"Component denominator: {denominator}")

    analytical_C = (p1 * p2 * p3) / denominator
    print(f"\nValue of the limit constant C = {analytical_C}\n")

    # Part 2: Numerically verify the integral to confirm the analytical calculation.
    def sech(x):
        return 1.0 / np.cosh(x)

    def integrand(z):
        # This is the function whose definite integral from -inf to inf gives C.
        return (np.pi / 2.0) * (sech(z / np.sqrt(2.0))**4)

    # Use scipy.integrate.quad for numerical integration.
    numerical_C, error = quad(integrand, -np.inf, np.inf)

    print("--- Numerical Verification ---")
    print(f"The numerically computed value of the integral is: {numerical_C}")
    print(f"Estimated error of the numerical integration is: {error}\n")

    # The fact that C is a finite positive number confirms that a=2 is attainable.
    # Since we established that a <= 2 for any solution, the largest possible value is 2.
    final_a = 2
    print(f"--- Conclusion ---")
    print(f"The analysis shows that a <= 2 for any solution, and a = 2 is attainable.")
    print(f"Therefore, the largest possible value for 'a' is: {final_a}")


solve()