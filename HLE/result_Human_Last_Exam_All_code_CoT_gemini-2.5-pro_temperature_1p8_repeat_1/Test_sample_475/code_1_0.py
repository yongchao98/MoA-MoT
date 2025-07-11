import numpy as np
from scipy.special import lambertw

def calculate_total_charge():
    """
    Calculates the total charge on the droplet based on the derived formula.

    This calculation proceeds under the assumption that the parameter q_i is 0,
    which makes the problem solvable and independent of the unspecified
    perturbation parameters epsilon, n, and m.
    """
    # Given constants
    sigma_0 = 7.43e-7  # units: e/nm
    R_0 = 30.0         # units: nm

    # Mathematical constants
    pi = np.pi
    
    # Calculate omega, which is the Lambert W function evaluated at 1.
    # We take the real part, as lambertw can return a complex number.
    omega = lambertw(1).real

    # Calculate the constant factor derived from the Lambert W function expression
    w_factor = omega / (1 + omega)**3

    # Final formula for the total charge Q: Q = pi^3 * sigma_0 * R_0 * [w / (1+w)^3]
    total_charge = pi**3 * sigma_0 * R_0 * w_factor

    # Output the steps of the final calculation with numerical values
    print("Final equation for the total charge Q:")
    print(f"Q = π^3 * σ₀ * R₀ * [ω / (1 + ω)^3]")
    print("\nSubstituting the numerical values:")
    print(f"π = {pi:.5f}")
    print(f"σ₀ = {sigma_0:.2e} e/nm")
    print(f"R₀ = {R_0:.1f} nm")
    print(f"ω = W(1) = {omega:.5f}")
    print("\nCalculation:")
    print(f"Q = ({pi:.5f})^3 * ({sigma_0:.2e}) * ({R_0:.1f}) * [{omega:.5f} / (1 + {omega:.5f})^3]")
    print(f"Q = ({pi**3:.5f}) * ({sigma_0:.2e}) * ({R_0:.1f}) * [{w_factor:.5f}]")
    
    # Print the final result
    print(f"\nThe calculated total charge on the droplet is:")
    print(f"Q = {total_charge:.4e} e")

if __name__ == "__main__":
    calculate_total_charge()
