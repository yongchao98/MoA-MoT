import numpy as np
from scipy.integrate import quad

def solve_and_explain():
    """
    This script solves the problem by analyzing a specific solution and using numerical
    integration to verify its properties, then generalizing the result.
    """

    # --- Step 1: Define the problem and a specific solution ---
    print("Step 1: Analyze the PDE and find a suitable solution.")
    print("The PDE is Delta u = u^3 - u, derived from the potential W(t) = 1/4 * (1-t^2)^2.")
    print("A one-dimensional solution that satisfies the condition |u|<1 is u(x,y,z) = tanh(x/sqrt(2)).")
    print("For this solution, the squared gradient is |nabla u|^2 = 1/2 * sech^4(x/sqrt(2)), where sech(z) = 1/cosh(z).")
    print("-" * 60)

    # --- Step 2: Analyze the energy integral for this solution ---
    print("Step 2: Calculate the growth of the energy integral E(R) = integral_{B_R} |nabla u|^2 dV.")
    print("The integral can be computed as:")
    print("E(R) = integral from x=-R to R of [ Area(disk_x) * |nabla u(x)|^2 ] dx")
    print("E(R) = integral from x=-R to R of [ pi * (R^2 - x^2) * 1/2 * sech^4(x/sqrt(2)) ] dx")
    print("\nWe will now compute this integral numerically for increasing R and check the ratio E(R) / R^2.")
    print("-" * 60)

    def sech(z):
        return 1 / np.cosh(z)

    def integrand(x, R):
        # Integrand for the energy calculation.
        # For a fixed x, the integration over y and z is over a disk of radius sqrt(R^2 - x^2).
        # The area of this disk is pi * (R^2 - x^2).
        nabla_u_sq = 0.5 * sech(x / np.sqrt(2))**4
        area_disk = np.pi * (R**2 - x**2)
        return nabla_u_sq * area_disk

    def calculate_energy(R):
        # Integrate the function from -R to R.
        result, _ = quad(integrand, -R, R, args=(R,))
        return result

    # --- Step 3: Perform numerical verification ---
    print("Step 3: Numerical verification of the energy growth.")
    radii = [10.0, 20.0, 50.0, 100.0]
    print(f"{'R':<10}{'Energy E(R)':<25}{'E(R) / R^2':<20}")
    print("-" * 60)
    for R in radii:
        energy = calculate_energy(R)
        ratio = energy / (R**2)
        print(f"{R:<10.1f}{energy:<25.4f}{ratio:<20.4f}")

    analytical_constant = (2 * np.sqrt(2) * np.pi) / 3
    print(f"\nThe analytical limit for E(R)/R^2 is (2*sqrt(2)*pi)/3 ≈ {analytical_constant:.4f}.")
    print("The numerical results confirm that for this solution, the energy E(R) grows proportionally to R^2.")
    print("-" * 60)

    # --- Step 4: Conclude the final answer ---
    print("Step 4: Determine the largest possible value of 'a'.")
    print("For our specific solution, E(R) ≈ C * R^2. Therefore, the condition")
    print("  liminf_{R->inf} R^{-a} E(R) > 0")
    print("is satisfied for a = 2. This shows that the largest possible 'a' is at least 2.")
    
    print("\nGeneral theory for the Allen-Cahn equation suggests that the energy of any valid solution")
    print("is concentrated on (n-1)-dimensional surfaces. In R^3, this means 2D surfaces.")
    print("The area of such a surface inside a ball B_R grows at most like R^2.")
    print("This implies that for ANY solution, E(R) <= C * R^2 for some constant C.")
    print("If a > 2, then R^{-a} E(R) <= C * R^{2-a}, which goes to 0 as R -> inf.")
    print("Therefore, the liminf cannot be positive if a > 2.")
    
    print("\nCombining these two points, the largest possible value for 'a' must be 2.")
    print("-" * 60)
    
    a = 2
    print(f"The final answer for the largest possible value of a is: {a}")

if __name__ == '__main__':
    solve_and_explain()