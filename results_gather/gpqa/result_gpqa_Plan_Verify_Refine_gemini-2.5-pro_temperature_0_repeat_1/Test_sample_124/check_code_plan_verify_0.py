import math

def check_answer():
    """
    This function checks the correctness of the provided answer for the energy and degeneracy
    of the third excited state of a 3D isotropic harmonic oscillator.
    """
    
    # The question describes a 3D isotropic harmonic oscillator.
    # The energy eigenvalues for this system are given by the formula:
    # E_n = (n + 3/2) * hbar * omega
    # where n = 0, 1, 2, ... is the principal quantum number.
    
    # The degeneracy (number of linearly independent eigenfunctions) for an energy level E_n is:
    # g_n = (n + 2) * (n + 1) / 2
    
    # The question asks for the third excited state.
    # Ground state: n = 0
    # First excited state: n = 1
    # Second excited state: n = 2
    # Third excited state: n = 3
    n = 3
    
    # The provided answer (Option A) states:
    # Energy = (9/2) * hbar * omega
    # Degeneracy = 10
    
    # Let's extract the expected values from the answer.
    # The energy is given as a factor of hbar * omega.
    expected_energy_factor = 9.0 / 2.0
    expected_degeneracy = 10
    
    # --- Perform Calculations based on the formulas ---
    
    # 1. Calculate the energy factor for n=3.
    calculated_energy_factor = n + 1.5
    
    # 2. Calculate the degeneracy for n=3.
    # The result must be an integer, so we use integer division.
    calculated_degeneracy = (n + 2) * (n + 1) // 2
    
    # --- Verify the calculated values against the answer ---
    
    energy_matches = math.isclose(calculated_energy_factor, expected_energy_factor)
    degeneracy_matches = (calculated_degeneracy == expected_degeneracy)
    
    if energy_matches and degeneracy_matches:
        return "Correct"
    else:
        error_messages = []
        if not energy_matches:
            error_messages.append(
                f"Constraint violated: Energy value. "
                f"For the third excited state (n=3), the calculated energy is ({n} + 3/2)hbar*omega = {calculated_energy_factor}*hbar*omega. "
                f"The answer provides an energy factor of {expected_energy_factor}, which is inconsistent."
            )
        if not degeneracy_matches:
            error_messages.append(
                f"Constraint violated: Degeneracy. "
                f"For n=3, the calculated degeneracy is (({n}+2)*({n}+1))/2 = {calculated_degeneracy}. "
                f"The answer provides a degeneracy of {expected_degeneracy}, which is inconsistent."
            )
        return "\n".join(error_messages)

# Execute the check and print the result.
# This is the code block that will be run to check the answer.
result = check_answer()
print(result)