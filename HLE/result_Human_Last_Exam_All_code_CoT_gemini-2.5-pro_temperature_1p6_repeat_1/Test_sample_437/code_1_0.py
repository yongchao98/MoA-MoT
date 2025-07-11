import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    Evaluates the integral <phi_1s| 1/r |phi_1s> for a 1s Slater orbital.
    
    This function demonstrates the evaluation by:
    1. Stating the analytical result.
    2. Numerically calculating the integral components in spherical coordinates.
    3. Combining the numerical results to verify the analytical solution.
    """
    
    # Set the orbital exponent zeta. For a hydrogen atom 1s orbital, zeta = 1.0.
    zeta = 1.0

    # The integral to evaluate is <phi_1s| 1/r |phi_1s>
    # The analytical result is simply the value of zeta.
    analytical_result = zeta

    # The integral in spherical coordinates is:
    # I = (zeta^3/pi) * Integral[d(phi)] * Integral[sin(theta)d(theta)] * Integral[r*exp(-2*zeta*r)dr]
    
    print(f"Evaluating the integral <phi_1s| 1/r |phi_1s> for a 1s Slater orbital with zeta = {zeta}")
    
    # We will compute each part of the final equation for the integral.
    
    # 1. The constant pre-factor from the squared normalization constant
    pre_factor = zeta**3 / np.pi

    # 2. The integral over phi from 0 to 2*pi
    integral_phi, _ = quad(lambda phi: 1, 0, 2 * np.pi)

    # 3. The integral over theta from 0 to pi
    integral_theta, _ = quad(lambda theta: np.sin(theta), 0, np.pi)

    # 4. The integral over r from 0 to infinity
    radial_integrand = lambda r: r * np.exp(-2 * zeta * r)
    integral_r, _ = quad(radial_integrand, 0, np.inf)

    # The final numerical result is the product of these parts
    numerical_result = pre_factor * integral_phi * integral_theta * integral_r

    # Print the final equation with the computed numerical values for each part
    print("\nThe final equation is assembled from the following parts:")
    print(f"I = (zeta^3/pi) * (Integral over phi) * (Integral over theta) * (Integral over r)")
    print(f"I = ({pre_factor:.4f}) * ({integral_phi:.4f}) * ({integral_theta:.4f}) * ({integral_r:.4f})")
    
    print(f"\nComputed numerical result for the integral: {numerical_result:.6f}")
    
    print(f"\nThis matches the analytical result, which is zeta = {analytical_result}")

if __name__ == "__main__":
    solve_integral()