import scipy.constants as const
import numpy as np

def check_answer():
    """
    Calculates the minimum potential energy of the system and checks it against the provided answer.
    """
    # --- Define Constants and Parameters ---
    # Physical constants from scipy for high precision
    k = const.k  # Coulomb's constant in N m^2 / C^2
    e = const.e  # Elementary charge in C

    # Parameters from the problem statement
    num_sphere_charges = 12
    charge_q = 2 * e
    radius_r = 2.0  # meters

    # The dimensionless energy constant for the icosahedral arrangement (N=12 Thomson problem).
    # This value is well-established in the literature.
    E_12 = 49.165

    # --- Calculation ---
    # The total potential energy is the sum of two components:
    # 1. U_cs: Interaction between the central charge and the 12 sphere charges.
    # 2. U_ss: Mutual interaction among the 12 sphere charges in their minimum energy configuration (icosahedron).
    # U_total = U_cs + U_ss
    # U_total = (12 * k * q^2 / r) + (E_12 * k * q^2 / r)
    # U_total = (12 + E_12) * (k * q^2 / r)
    
    total_energy_factor = num_sphere_charges + E_12
    base_energy_unit = (k * charge_q**2) / radius_r
    
    calculated_energy = total_energy_factor * base_energy_unit

    # --- Verification ---
    # The provided answer is C, which corresponds to the value 2.822 x 10^-26 J.
    expected_energy = 2.822e-26

    # Check if the calculated value matches the expected value.
    # A relative tolerance of 1e-3 is appropriate, as the options are given to 3 decimal places
    # on the mantissa (e.g., 2.822), which implies a precision of 1 part in 1000.
    if np.isclose(calculated_energy, expected_energy, rtol=1e-3):
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {expected_energy:.4e} J. "
                f"However, the calculated minimum energy is {calculated_energy:.4e} J. "
                f"The calculation is based on the formula U_total = (12 + E_12) * (k * q^2 / r), "
                f"with E_12 = 49.165 for the icosahedron configuration, which is the correct physical model.")

# Run the check
result = check_answer()
print(result)