import math

def check_diatomic_molecule_momentum():
    """
    Calculates the momentum of a photon absorbed by a diatomic molecule
    undergoing a rovibrational transition from the ground state.
    """
    # --- Constants ---
    h_bar = 1.054571817e-34  # Reduced Planck constant (J*s)
    c = 2.99792458e8         # Speed of light (m/s)
    amu_to_kg = 1.660539e-27   # Atomic mass unit to kg conversion
    angstrom_to_m = 1e-10      # Angstrom to meter conversion

    # --- Given parameters from the question ---
    Mx_amu = 20.0
    My_amu = 2.0
    R_angstrom = 2.0
    omega = 4.0e14  # rad/s

    # --- The answer to check (Option C) ---
    # The provided final answer is <<<C>>>, which corresponds to p = 1.4*10^(-28) N*s
    expected_p = 1.4e-28

    # --- Step-by-step calculation ---

    # 1. Convert inputs to SI units
    Mx_kg = Mx_amu * amu_to_kg
    My_kg = My_amu * amu_to_kg
    R_m = R_angstrom * angstrom_to_m

    # 2. Calculate reduced mass (mu)
    # Formula: mu = (m1 * m2) / (m1 + m2)
    mu = (Mx_kg * My_kg) / (Mx_kg + My_kg)

    # 3. Calculate moment of inertia (I)
    # Formula: I = mu * R^2
    I = mu * (R_m ** 2)

    # 4. Calculate the energy difference (delta_E)
    # The transition is from the ground state (v=0, J=0) to the next allowed
    # rovibrational state (v=1, J=1).
    # E_initial = (1/2) * h_bar * omega
    # E_final = (3/2) * h_bar * omega + h_bar^2 / I
    # delta_E = E_final - E_initial = h_bar * omega + h_bar^2 / I
    vibrational_energy_term = h_bar * omega
    rotational_energy_term = (h_bar ** 2) / I
    delta_E = vibrational_energy_term + rotational_energy_term

    # 5. Calculate the photon's momentum (p)
    # Formula: p = delta_E / c
    calculated_p = delta_E / c

    # --- Verification ---
    # Check if the calculated momentum is close to the expected value from option C.
    # A relative tolerance of 2% is used to account for rounding of constants.
    tolerance = 0.02
    if abs(calculated_p - expected_p) / expected_p < tolerance:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer C (p = {expected_p:.2e} N*s) is not correct.\n"
                f"The calculated momentum is p = {calculated_p:.3e} N*s.\n"
                f"Details:\n"
                f"  - Reduced mass (μ): {mu:.4e} kg\n"
                f"  - Moment of inertia (I): {I:.4e} kg·m²\n"
                f"  - Energy difference (ΔE): {delta_E:.4e} J\n")

# Execute the check and print the result
result = check_diatomic_molecule_momentum()
print(result)