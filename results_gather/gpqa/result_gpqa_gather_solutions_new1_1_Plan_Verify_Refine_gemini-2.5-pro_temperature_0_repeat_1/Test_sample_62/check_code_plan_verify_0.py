import numpy as np
from scipy import constants

def check_diatomic_molecule_momentum():
    """
    Calculates the momentum of a photon absorbed by a diatomic molecule
    for a specific rovibrational transition and checks it against the given answer.
    """
    # --- 1. Define Physical Constants and Given Parameters ---
    # Use high-precision values from scipy.constants
    h_bar = constants.hbar  # Reduced Planck constant (J·s)
    c = constants.c         # Speed of light (m/s)
    amu_to_kg = constants.physical_constants['atomic mass constant'][0]

    # Given parameters from the question
    Mx_amu = 20.0
    My_amu = 2.0
    R_angstrom = 2.0
    omega = 4.0e14  # rad/s

    # --- 2. Convert all parameters to base SI units ---
    Mx_kg = Mx_amu * amu_to_kg
    My_kg = My_amu * amu_to_kg
    R_m = R_angstrom * 1e-10

    # --- 3. Identify the transition and calculate the energy difference (ΔE) ---
    # The molecule starts in the fundamental (ground) state: v=0, J=0.
    # The selection rules for photon absorption are Δv=+1 and ΔJ=±1.
    # From J=0, only ΔJ=+1 is possible.
    # Therefore, the final state is v=1, J=1.
    
    # The energy of a state is E(v, J) = ħω(v + 1/2) + (ħ²/2I) * J(J + 1).
    # The energy difference ΔE = E(1,1) - E(0,0) simplifies to:
    # ΔE = [ħω(3/2) + (ħ²/I)] - [ħω(1/2)] = ħω + ħ²/I

    # To calculate ΔE, we first need the moment of inertia (I).
    
    # Calculate the reduced mass (μ)
    mu = (Mx_kg * My_kg) / (Mx_kg + My_kg)

    # Calculate the moment of inertia (I)
    I = mu * (R_m ** 2)

    # Calculate the total transition energy (ΔE)
    delta_E = (h_bar * omega) + ((h_bar ** 2) / I)

    # --- 4. Calculate the photon's momentum (p) ---
    p_calculated = delta_E / c

    # --- 5. Compare the calculated value with the provided answer ---
    # The final answer from the analysis is <<<B>>>, which corresponds to p = 1.4 * 10^(-28) N*s.
    answer_value = 1.4e-28

    # Check if the calculated value is close to the answer's value.
    # A relative tolerance of 2% is reasonable given the rounding in the options.
    if np.isclose(p_calculated, answer_value, rtol=0.02):
        return "Correct"
    else:
        # If incorrect, provide a detailed reason.
        reason = (
            f"Incorrect. The provided answer corresponds to a momentum of {answer_value:.2e} N*s.\n"
            f"The calculated momentum is {p_calculated:.4e} N*s.\n\n"
            "Calculation Breakdown:\n"
            f"1. Reduced Mass (μ): {mu:.4e} kg\n"
            f"2. Moment of Inertia (I): {I:.4e} kg·m²\n"
            f"3. Transition Energy (ΔE): {delta_E:.4e} J\n"
            f"4. Calculated Momentum (p = ΔE/c): {p_calculated:.4e} N*s"
        )
        return reason

# Execute the check and print the result.
result = check_diatomic_molecule_momentum()
print(result)