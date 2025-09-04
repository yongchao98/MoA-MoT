import re

def check_quantum_harmonic_oscillator_answer():
    """
    Checks the correctness of the answer for the 3D isotropic harmonic oscillator problem.
    """
    # --- Define Problem Constraints ---
    # The question asks for the "third excited state".
    # Ground state is n=0, first excited is n=1, second is n=2, third is n=3.
    principal_quantum_number = 3

    # --- Calculate Correct Theoretical Values ---
    # Energy formula: E_n = (n + 3/2) * hbar * omega
    # We check the numerical factor in front of hbar*omega.
    correct_energy_factor = principal_quantum_number + 1.5  # Should be 4.5

    # Degeneracy formula: g_n = (n+1)(n+2)/2
    correct_degeneracy = (principal_quantum_number + 1) * (principal_quantum_number + 2) // 2 # Should be 10

    # --- Parse the LLM's Answer ---
    llm_answer_text = "<<<D>>>"
    
    # Define the options as given in the question
    options = {
        'A': {'energy_str': "(9/2) \\hbar \\omega", 'degeneracy': 3},
        'B': {'energy_str': "11 \\pi^2 \\hbar^2 / (2m r^2)", 'degeneracy': 10},
        'C': {'energy_str': "11 \\pi^2 \\hbar^2 / (2m r^2)", 'degeneracy': 3},
        'D': {'energy_str': "(9/2) \\hbar \\omega", 'degeneracy': 10}
    }

    # Extract the chosen option letter
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Failure: Could not parse the chosen option from the answer string."
    chosen_option_letter = match.group(1)
    
    chosen_option = options[chosen_option_letter]

    # --- Verify the Chosen Option ---
    # Check 1: Verify the energy formula type and value.
    # The correct energy is proportional to hbar*omega, not 1/r^2.
    if "hbar" in chosen_option['energy_str'] and "omega" in chosen_option['energy_str']:
        # This is the correct formula type. Now check the factor.
        # For (9/2), the factor is 4.5
        if chosen_option['energy_str'].startswith("(9/2)"):
            chosen_energy_factor = 9.0 / 2.0
            if chosen_energy_factor != correct_energy_factor:
                return (f"Incorrect. The energy value is wrong. "
                        f"For n={principal_quantum_number}, the energy should be ({int(2*correct_energy_factor)}/2)hbar*omega, "
                        f"but option {chosen_option_letter} has a factor of {chosen_energy_factor}.")
        else:
            # This case handles if there was another option with the correct formula type but wrong value.
            return f"Incorrect. The energy value in option {chosen_option_letter} is not correct."
    else:
        # The energy formula is for a different system (e.g., particle in a spherical well)
        return (f"Incorrect. The energy formula in option {chosen_option_letter} is for the wrong physical system. "
                f"The potential V(r) ~ r^2 corresponds to a harmonic oscillator.")

    # Check 2: Verify the degeneracy.
    if chosen_option['degeneracy'] != correct_degeneracy:
        return (f"Incorrect. The degeneracy is wrong. "
                f"For n={principal_quantum_number}, the degeneracy should be {correct_degeneracy}, "
                f"but option {chosen_option_letter} has a degeneracy of {chosen_option['degeneracy']}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_quantum_harmonic_oscillator_answer()
print(result)