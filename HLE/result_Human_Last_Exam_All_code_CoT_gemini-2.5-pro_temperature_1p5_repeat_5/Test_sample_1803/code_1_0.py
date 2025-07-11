import numpy as np
from scipy.constants import hbar, e, epsilon_0, m_e, pi

def calculate_ground_state_energy_shift(m, omega_0, R):
    """
    Calculates the leading term of the ground state energy shift for two
    interacting quantum harmonic oscillators.

    The formula for the energy shift (Delta E_0) is:
    Delta E_0 = - (3 * hbar * e**4) / (64 * pi**2 * epsilon_0**2 * m**2 * omega_0**3 * R**6)

    Args:
        m (float): Mass of the oscillator in kg.
        omega_0 (float): Angular frequency of the oscillator in rad/s.
        R (float): Distance between the oscillators in meters.

    Returns:
        float: The ground state energy shift in Joules.
    """
    # Calculate the numerator of the formula
    numerator = 3 * hbar * e**4

    # Calculate the denominator of the formula
    denominator = 64 * pi**2 * epsilon_0**2 * m**2 * omega_0**3 * R**6
    
    # Calculate the total energy shift
    energy_shift = -numerator / denominator

    # Print the equation with the calculated numbers
    print("Calculating the ground state energy shift Delta E_0:")
    print(f"Formula: Delta E_0 = - (3 * hbar * e^4) / (64 * pi^2 * epsilon_0^2 * m^2 * omega_0^3 * R^6)")
    print("\n--- Numerical Values ---")
    print(f"Numerator (3 * hbar * e^4): {numerator:.4e}")
    print(f"Denominator (64 * pi^2 * epsilon_0^2 * m^2 * omega_0^3 * R^6): {denominator:.4e}")
    print("------------------------")
    print(f"Final Energy Shift (Delta E_0) = -({numerator:.4e}) / ({denominator:.4e})")
    print(f"Delta E_0 = {energy_shift:.4e} Joules")
    # Also print in electron-volts (eV) for better intuition
    energy_shift_eV = energy_shift / e
    print(f"Delta E_0 = {energy_shift_eV:.4e} eV")
    
    return energy_shift

if __name__ == '__main__':
    # --- Example Parameters ---
    # We can use properties of a hydrogen atom for a reasonable scale.
    # Mass (electron mass)
    mass = m_e
    # A typical atomic transition energy is a few eV. Let's use E = 10 eV.
    # E = hbar * omega_0 => omega_0 = E / hbar
    energy_eV = 10
    angular_frequency = (energy_eV * e) / hbar
    # Distance between oscillators (e.g., 1 nanometer)
    distance = 1e-9  # 1 nm

    print("--- Input Parameters ---")
    print(f"Mass (m): {mass:.4e} kg")
    print(f"Angular Frequency (omega_0): {angular_frequency:.4e} rad/s")
    print(f"Distance (R): {distance:.4e} m\n")

    calculate_ground_state_energy_shift(m=mass, omega_0=angular_frequency, R=distance)