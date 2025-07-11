import numpy as np
from scipy.integrate import quad

def solve():
    """
    This script analyzes the growth rate of the integral of |nabla u|^2
    for the specific 1D solution u(x, y, z) = tanh(x/sqrt(2)) to the PDE
    Delta u = u^3 - u.
    
    The goal is to find the largest 'a' for which the following holds:
    liminf_{R->inf} R^(-a) * Integral_{B_R} |nabla u|^2 > 0
    """

    # For u(x) = tanh(x/sqrt(2)), the squared gradient is |nabla u|^2 = 1/2 * sech^4(x/sqrt(2)).
    # The integral over a ball B_R of radius R is given by:
    # I(R) = Integral from -R to R of [ pi * (R^2 - x^2) * 1/2 * sech^4(x/sqrt(2)) ] dx
    # because the cross-sectional area at x is pi*(R^2-x^2).
    
    def integrand(x, R):
        # Using the definition of sech(v) = 1/cosh(v)
        sech_val = 1 / np.cosh(x / np.sqrt(2))
        return np.pi * (R**2 - x**2) * 0.5 * sech_val**4

    # We test the scaling for large R. We expect I(R) to scale with R^2.
    # So, the ratio I(R) / R^2 should converge to a constant C.
    R_values = [10, 20, 50, 100]

    print("Analyzing the integral I(R) = Integral_{B_R} |nabla u|^2 dV for u(x,y,z)=tanh(x/sqrt(2)).")
    print("We compute the ratio I(R) / R^2 to find the scaling behavior.\n")

    for R in R_values:
        integral_val, _ = quad(integrand, -R, R, args=(R,))
        ratio = integral_val / R**2
        print(f"For R = {R}:")
        print(f"  Integral I(R) = {integral_val:.4f}")
        print(f"  Ratio I(R) / R^2 = {ratio:.4f}")

    # The analytical limit C = lim_{R->inf} I(R) / R^2 is (4 * pi) / (3 * sqrt(2)).
    C_analytical = (4 * np.pi) / (3 * np.sqrt(2))
    print(f"\nAs R increases, the ratio converges to C, analytically calculated as {C_analytical:.4f}.")
    print("This shows that for this solution, the integral grows proportionally to R^2.")
    
    print("\nThe problem is to find the largest 'a' satisfying the condition.")
    print("Substituting the growth rate, we get:")
    print("liminf_{R->inf} R^(-a) * (C * R^2) > 0")
    print("This simplifies to liminf_{R->inf} C * R^(2-a) > 0.")
    
    # Final step: determine a
    final_exponent = 2
    print(f"\nFor the limit to be positive, the exponent of R must be non-negative:")
    print(f"{final_exponent} - a >= 0  =>  a <= {final_exponent}")
    
    print("\nThis specific solution demonstrates that a=2 is possible.")
    print("Theory suggests that this is the maximal growth rate for any such solution.")
    print("Therefore, the largest possible value for 'a' is 2.")

solve()
<<<2>>>