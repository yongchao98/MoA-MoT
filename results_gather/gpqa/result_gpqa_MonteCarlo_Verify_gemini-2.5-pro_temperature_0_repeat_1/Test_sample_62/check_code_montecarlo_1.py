import math
from scipy import constants

def check_diatomic_molecule_transition():
    """
    This function calculates the momentum of a photon required to excite a diatomic
    molecule from its ground state to the next lowest possible energy state,
    and checks it against the provided answer.
    """
    # --- Given Parameters from the Question ---
    Mx_amu = 20.0  # Mass of atom X in amu
    My_amu = 2.0   # Mass of atom Y in amu
    R_angstrom = 2.0  # Molecular bond length in angstroms
    omega = 4.0 * 10**14  # Angular frequency of vibration in rad/s

    # --- Provided Answer ---
    # The question's answer is D, with a value of 1.4e-28 N*s
    answer_value = 1.4e-28

    # --- Physical Constants from scipy.constants for high precision ---
    amu_to_kg = constants.atomic_mass  # Conversion factor from amu to kg
    hbar = constants.hbar              # Reduced Planck constant in J*s
    c = constants.c                    # Speed of light in m/s

    # --- Step 1: Convert all inputs to SI units ---
    Mx_kg = Mx_amu * amu_to_kg
    My_kg = My_amu * amu_to_kg
    R_m = R_angstrom * 1e-10

    # --- Step 2: Calculate intermediate physical properties ---

    # Reduced mass (μ)
    try:
        mu = (Mx_kg * My_kg) / (Mx_kg + My_kg)
    except ZeroDivisionError:
        return "Error: Sum of masses cannot be zero."

    # Moment of inertia (I)
    I = mu * R_m**2

    # --- Step 3: Calculate the energy of the transition (ΔE) ---
    # The transition is from the ground state (v=0, J=0) to the next lowest
    # accessible state via photon absorption. Selection rules (Δv=+1, ΔJ=+1)
    # dictate the final state is (v=1, J=1).
    # ΔE = E(1,1) - E(0,0) = (ħω) + (ħ²/I)
    
    energy_vibrational = hbar * omega
    energy_rotational = (hbar**2) / I
    delta_E = energy_vibrational + energy_rotational

    # --- Step 4: Calculate the photon's momentum (p) ---
    calculated_p = delta_E / c

    # --- Step 5: Check correctness of the answer ---
    # We use math.isclose to check if the calculated value is close to the
    # answer's value, allowing for minor rounding differences. A relative
    # tolerance of 2% is reasonable for this type of problem.
    if math.isclose(calculated_p, answer_value, rel_tol=0.02):
        return "Correct"
    else:
        return (f"Incorrect. The calculated momentum is {calculated_p:.4e} N*s. "
                f"The provided answer D corresponds to a value of {answer_value:.4e} N*s. "
                f"The calculated value differs by {abs(calculated_p - answer_value)/answer_value:.2%} which is outside the 2% tolerance.")

# Execute the check and print the result
result = check_diatomic_molecule_transition()
print(result)