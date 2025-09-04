import scipy.constants as const
import numpy as np

# This script verifies the correct answer for the physics problem by calculating it from first principles.
# It checks if the calculation matches the expected correct option derived from the consensus of the provided answers.

def check_final_answer():
    """
    Calculates the minimum potential energy and verifies the correctness of the consensus answer.
    Returns "Correct" if the calculation matches the expected result, otherwise returns a reason for the discrepancy.
    """
    try:
        # --- 1. Define Constants and Parameters ---
        # Use precise physical constants from the scipy library.
        k = const.physical_constants['Coulomb constant'][0]  # Coulomb's constant (N·m²/C²)
        e = const.e                                          # Elementary charge (C)

        # Parameters from the problem statement.
        q = 2 * e      # Charge of each particle
        r = 2.0        # Radius of the sphere (m)
        N_sphere = 12  # Number of charges on the sphere

        # For the Thomson problem with N=12, the minimum energy configuration is an icosahedron.
        # The dimensionless energy constant E_12 represents the sum of inverse chord lengths
        # for a unit-radius icosahedron. Its value is well-established as approximately 49.165.
        E_12 = 49.165

        # --- 2. Calculate the Total Minimum Potential Energy ---
        # The total energy U_total is the sum of the interaction between the center and sphere charges (U_cs)
        # and the mutual interaction between the sphere charges (U_ss).
        # U_total = U_cs + U_ss = (N_sphere * k * q**2 / r) + (E_12 * k * q**2 / r)
        # U_total = (N_sphere + E_12) * (k * q**2 / r)
        total_energy = (N_sphere + E_12) * (k * q**2 / r)

        # --- 3. Check Against the Provided Options ---
        # The options from the original question prompt are:
        # A) 5.645 x 10^-26
        # B) 7.056 x 10^-27
        # C) 2.822 x 10^-26
        # D) 122.330
        options = {
            'A': 5.645e-26,
            'B': 7.056e-27,
            'C': 2.822e-26,
            'D': 122.330
        }

        # Find which option key corresponds to the calculated value using a relative tolerance.
        correct_option = None
        for key, value in options.items():
            if np.isclose(total_energy, value, rtol=1e-3):
                correct_option = key
                break
        
        # --- 4. Final Verification ---
        # The consensus among the correct candidate answers and the final provided script is that
        # the result is ~2.822e-26 J, which corresponds to option C.
        # We check if our independent calculation confirms this.
        expected_option = 'C'

        if correct_option == expected_option:
            # The calculation confirms the consensus answer.
            return "Correct"
        elif correct_option is None:
            return (f"Incorrect. The calculated value {total_energy:.4e} J does not match any of the provided options "
                    f"within the tolerance. There may be an error in the question's options.")
        else:
            # The calculation does not confirm the consensus answer.
            return (f"Incorrect. The independent calculation yields a value of {total_energy:.4e} J, "
                    f"which corresponds to option '{correct_option}'. This does not match the expected "
                    f"correct option '{expected_option}'.")

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result.
print(check_final_answer())