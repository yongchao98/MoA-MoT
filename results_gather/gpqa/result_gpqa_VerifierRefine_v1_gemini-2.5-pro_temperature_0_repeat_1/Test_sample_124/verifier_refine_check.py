import math

def check_3d_harmonic_oscillator_answer():
    """
    This function checks the energy and degeneracy for the third excited state
    of a 3D isotropic harmonic oscillator.

    The energy levels are given by E_N = (N + 3/2) * hbar * omega, where N = 0, 1, 2, ...
    The degeneracy of level N is g_N = (N + 1) * (N + 2) / 2.

    The ground state is N=0.
    The first excited state is N=1.
    The second excited state is N=2.
    The third excited state is N=3.
    """
    
    # The question asks for the third excited state.
    # Ground state is n=0, 1st excited is n=1, 2nd is n=2, 3rd is n=3.
    principal_quantum_number_N = 3

    # --- 1. Calculate the theoretical energy ---
    # The energy is E = (N + 3/2) * hbar * omega.
    # We will check the numerical factor (N + 3/2).
    calculated_energy_factor = principal_quantum_number_N + 3/2

    # --- 2. Calculate the theoretical degeneracy ---
    # The degeneracy is g = (N + 1) * (N + 2) / 2.
    calculated_degeneracy = (principal_quantum_number_N + 1) * (principal_quantum_number_N + 2) / 2
    
    # Ensure the result is an integer
    calculated_degeneracy = int(calculated_degeneracy)

    # --- 3. Define the values from the provided answer (Option B) ---
    # Answer B states: Energy = (9/2) * hbar * omega, Degeneracy = 10
    expected_energy_factor = 9/2
    expected_degeneracy = 10

    # --- 4. Compare theoretical values with the answer's values ---
    energy_is_correct = (calculated_energy_factor == expected_energy_factor)
    degeneracy_is_correct = (calculated_degeneracy == expected_degeneracy)

    # --- 5. Formulate the final output ---
    if energy_is_correct and degeneracy_is_correct:
        return "Correct"
    else:
        error_messages = []
        if not energy_is_correct:
            error_messages.append(
                f"The energy value is incorrect. For the third excited state (N={principal_quantum_number_N}), "
                f"the energy should be ({principal_quantum_number_N} + 3/2) * hbar * omega = {calculated_energy_factor} * hbar * omega. "
                f"The answer provides an energy factor of {expected_energy_factor}."
            )
        if not degeneracy_is_correct:
            error_messages.append(
                f"The degeneracy is incorrect. For the third excited state (N={principal_quantum_number_N}), "
                f"the degeneracy is (N+1)(N+2)/2 = ({principal_quantum_number_N}+1)({principal_quantum_number_N}+2)/2 = {calculated_degeneracy}. "
                f"The answer provides a degeneracy of {expected_degeneracy}."
            )
        # Also check for the form of the energy, which is a common mistake.
        # The energy for a harmonic oscillator is proportional to hbar*omega, not hbar^2/(m*r^2).
        # Options A and C have the wrong form for the energy.
        if expected_energy_factor != 9/2:
             error_messages.append(
                "The functional form of the energy in options A and C, proportional to hbar^2/(m*r^2), "
                "is characteristic of a particle in an infinite spherical well, not a harmonic oscillator."
            )

        return "\n".join(error_messages)

# Execute the check and print the result
result = check_3d_harmonic_oscillator_answer()
print(result)