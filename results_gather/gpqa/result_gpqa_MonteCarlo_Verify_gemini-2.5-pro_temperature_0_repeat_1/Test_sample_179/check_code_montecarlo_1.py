import numpy as np

def check_solution():
    """
    Verifies the calculation for the minimum energy of the 13-particle system.
    The minimum energy configuration for the 12 shell particles is an icosahedron.
    This function calculates the total potential energy for this configuration and
    compares it to the provided answer.
    """
    # --- Define physical constants and system parameters from the question ---
    e = 1.602176634e-19  # Elementary charge in Coulombs
    k = 8.9875517923e9   # Coulomb's constant in N m^2 / C^2
    R = 2.0              # Radius of the sphere in meters
    N_shell = 12         # Number of particles on the shell
    q = 2 * e            # Charge of each particle

    # The value from the selected answer option D
    answer_value = 2.822e-26

    # --- Calculation ---
    # The total potential energy (U_total) is the sum of two components:
    # 1. U_center_shell: Interaction energy between the central charge and the 12 shell charges.
    # 2. U_shell_shell: Interaction energy among the 12 shell charges themselves.

    # 1. Calculate U_center_shell
    # This is constant as all 12 shell charges are at a fixed distance R from the center.
    U_center_shell = N_shell * k * q**2 / R

    # 2. Calculate U_shell_shell for the minimum energy configuration (icosahedron).
    # This can be written as U_shell_shell = (k * q^2 / R) * E_12, where E_12 is the
    # dimensionless energy constant for an icosahedron of 12 charges on a unit sphere.
    
    # The constant E_12 is derived from the geometry of the icosahedron.
    # It represents the sum of inverse distances for all N*(N-1)/2 = 66 pairs.
    phi = (1 + np.sqrt(5)) / 2  # The golden ratio
    
    # The sum of inverse distances from a single vertex to the other 11 on a unit sphere is:
    # (1/2) * [5 * phi * sqrt(1 + phi^2) + 1]
    # To get the total sum for all pairs, we multiply by N_shell and divide by 2 (to avoid double counting).
    # E_12 = (N_shell / 2) * (1/2) * [5 * phi * sqrt(1 + phi^2) + 1]
    E_12 = 6 * (1/2) * (5 * phi * np.sqrt(1 + phi**2) + 1)
    E_12 = 3 * (5 * phi * np.sqrt(1 + phi**2) + 1)
    
    U_shell_shell = (k * q**2 / R) * E_12

    # 3. Calculate the total minimum energy
    calculated_total_energy = U_center_shell + U_shell_shell

    # --- Verification ---
    # Check if the calculated energy matches the answer's value within a small tolerance
    # to account for rounding differences in physical constants.
    tolerance = 0.001  # 0.1% relative tolerance

    if abs(calculated_total_energy - answer_value) / calculated_total_energy < tolerance:
        return "Correct"
    else:
        return (f"Incorrect. The final answer is not consistent with the calculation.\n"
                f"The provided answer is {answer_value:.4e} J.\n"
                f"The calculated minimum energy is {calculated_total_energy:.4e} J.\n"
                f"Breakdown of calculated energy:\n"
                f"  - Central-Shell Energy: {U_center_shell:.4e} J\n"
                f"  - Shell-Shell Energy (Icosahedron): {U_shell_shell:.4e} J\n"
                f"The relative difference is {abs(calculated_total_energy - answer_value) / calculated_total_energy:.2%}, "
                f"which is outside the tolerance of {tolerance:.1%}.")

# Execute the check and print the result.
result = check_solution()
print(result)