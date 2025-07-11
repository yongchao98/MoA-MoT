import numpy as np
from scipy.optimize import root_scalar

def solve_energy_difference():
    """
    Calculates the energy difference between the first and second energy levels
    of a particle in an approximated 3D infinite spherical potential well.
    """
    # Physical constants
    m = 9.11e-31  # Mass of the particle in kg
    R = 3e-9      # Radius of the well in meters
    hbar = 1.05457e-34  # Reduced Planck constant in J·s
    eV_to_J = 1.60218e-19  # Conversion factor from eV to Joules

    # Find the zeros of the spherical Bessel functions
    # First zero of j_0(x) is z_10 = pi
    z10 = np.pi

    # First zero of j_1(x) is the first positive root of tan(x) = x
    # We find this root numerically. The root is in the interval (pi, 1.5*pi).
    def j1_zero_func(x):
        # np.tan(x) is large near 1.5*pi, but the root is well before it.
        # This check avoids overflow but isn't strictly necessary for the known root location.
        if np.cos(x) == 0:
            return np.inf
        return np.tan(x) - x

    # Search for the root in the interval (pi, 1.5*pi)
    result = root_scalar(j1_zero_func, bracket=[np.pi, 1.5 * np.pi - 0.01])
    z11 = result.root

    # The second zero of j_0(x) is z_20 = 2*pi.
    z20 = 2 * np.pi
    
    # The energy is proportional to z_nl^2.
    # z10^2 = pi^2 ~ 9.87
    # z11^2 ~ (4.493)^2 ~ 20.19
    # z20^2 = (2*pi)^2 ~ 39.48
    # The lowest two levels correspond to z10 and z11.

    # Calculate the pre-factor in the energy formula, converted to eV
    energy_factor_J = hbar**2 / (2 * m * R**2)
    energy_factor_eV = energy_factor_J / eV_to_J

    # Calculate the energy difference Delta E = E_2 - E_1
    delta_E_eV = energy_factor_eV * (z11**2 - z10**2)

    # Print the explanation and the final result
    print("This solution models the system as a particle in a 3D infinite spherical potential well.")
    print("The energy difference is calculated as Delta_E = (hbar^2 / (2*m*R^2)) * (z_11^2 - z_10^2)\n")
    print("Values used:")
    print(f"hbar = {hbar:.4e} J·s")
    print(f"m = {m:.2e} kg")
    print(f"R = {R:.1e} m")
    print(f"z_10 (first zero of j_0) = {z10:.4f}")
    print(f"z_11 (first zero of j_1) = {z11:.4f}\n")
    
    print("Step-by-step calculation of Delta_E:")
    print(f"Energy pre-factor (hbar^2 / (2*m*R^2)) = {energy_factor_eV:.6f} eV")
    print(f"Delta_E = {energy_factor_eV:.6f} eV * ({z11:.4f}^2 - {z10:.4f}^2)")
    print(f"Delta_E = {energy_factor_eV:.6f} eV * ({z11**2:.4f} - {z10**2:.4f})")
    print(f"Delta_E = {energy_factor_eV:.6f} eV * ({z11**2 - z10**2:.4f})")
    print(f"\nFinal energy difference Delta_E = {delta_E_eV:.4f} eV")
    
    # Final answer in specified format
    print(f"\n<<<{delta_E_eV:.4f}>>>")

solve_energy_difference()