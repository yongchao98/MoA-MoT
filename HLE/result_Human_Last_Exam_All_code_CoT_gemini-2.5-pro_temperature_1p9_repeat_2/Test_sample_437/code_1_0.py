import numpy as np
from scipy.integrate import quad

def evaluate_integral(zeta):
    """
    Evaluates the integral <phi_1s| 1/r |phi_1s> for a given zeta.

    The integral is broken down into three parts based on the analytical derivation:
    Result = (Normalization Constant)^2 * (Angular Integral) * (Radial Integral)

    Args:
        zeta (float): The orbital exponent of the 1s Slater orbital.
    
    Returns:
        float: The numerically calculated value of the integral.
    """

    # 1. Normalization constant squared for a 1s STO is N^2 = zeta^3 / pi
    norm_constant_squared = (zeta**3) / np.pi

    # 2. The angular part of the integral in spherical coordinates is 4*pi
    # Integral(d(phi)) from 0 to 2*pi is 2*pi
    # Integral(sin(theta)d(theta)) from 0 to pi is 2
    angular_integral = 4 * np.pi

    # 3. The radial part of the integral is Integral( r * exp(-2*zeta*r) dr ) from r=0 to infinity.
    # We define the integrand and compute it numerically using scipy.integrate.quad.
    radial_integrand = lambda r: r * np.exp(-2 * zeta * r)
    radial_integral_val, _ = quad(radial_integrand, 0, np.inf)

    # Combine the parts to get the final result
    final_result = norm_constant_squared * angular_integral * radial_integral_val
    
    # --- Output the step-by-step evaluation ---
    print(f"Evaluation of <phi_1s| 1/r |phi_1s> for a Slater orbital with zeta = {zeta}")
    print("-" * 70)
    print("The final equation is: Result = (N^2) * (Angular Integral) * (Radial Integral)")
    print(f"1. Normalization constant squared (N^2 = zeta^3/pi): {norm_constant_squared:.6f}")
    print(f"2. Angular Integral (4*pi): {angular_integral:.6f}")
    print(f"3. Radial Integral (from 0 to inf of r*exp(-2*zeta*r) dr): {radial_integral_val:.6f}")
    print("-" * 70)

    # Final calculation showing how the numerical components combine
    print(f"Final Numerical Result = {norm_constant_squared:.6f} * {angular_integral:.6f} * {radial_integral_val:.6f} = {final_result:.6f}")
    
    # Compare with the simple analytical result, which is just zeta
    print(f"\nThis numerical result matches the analytical result, which is exactly zeta = {zeta:.6f}")


# --- Main execution ---
if __name__ == "__main__":
    # You can change this value to test with other orbital exponents
    # For a Hydrogen atom in a minimal basis set (STO-1G), zeta is 1.0.
    # For a STO-3G basis set, the effective zeta for Hydrogen is ~1.24.
    zeta_value = 2.5
    evaluate_integral(zeta_value)
