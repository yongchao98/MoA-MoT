import math

def check_3d_harmonic_oscillator_answer():
    """
    This function checks the energy and degeneracy for the third excited state
    of a 3D isotropic quantum harmonic oscillator.

    The problem states:
    - System: Particle in a 3D isotropic potential V(r) = 1/2 * m * omega^2 * r^2.
    - Question: What is the energy and degeneracy of the third excited state?

    The provided answer to check is from option B:
    - Energy: (9/2) * hbar * omega
    - Degeneracy: 10
    """

    # The principal quantum number 'n' for the excited states is defined as:
    # n = 0: Ground state
    # n = 1: First excited state
    # n = 2: Second excited state
    # n = 3: Third excited state
    n = 3

    # --- Theoretical Calculation ---

    # 1. Calculate the expected energy.
    # The energy eigenvalues for a 3D isotropic harmonic oscillator are given by:
    # E_n = (n + 3/2) * hbar * omega
    # For the third excited state (n=3):
    # E_3 = (3 + 3/2) * hbar * omega = (9/2) * hbar * omega
    expected_energy_numerator = 2 * n + 3
    expected_energy_denominator = 2

    # 2. Calculate the expected degeneracy.
    # The degeneracy g_n of the energy level E_n is given by the formula:
    # g_n = (n + 1) * (n + 2) / 2
    # For the third excited state (n=3):
    # g_3 = (3 + 1) * (3 + 2) / 2 = 4 * 5 / 2 = 10
    expected_degeneracy = (n + 1) * (n + 2) // 2

    # --- Parse and Verify the Provided Answer (Option B) ---
    # Answer B provides: Energy = (9/2) hbar omega, Degeneracy = 10
    answer_energy_str = "(9/2) \\hbar \\omega"
    answer_degeneracy = 10

    # Verify the energy value by parsing the fraction
    try:
        # Extract the fraction string, e.g., "(9/2)"
        frac_str = answer_energy_str.split(' ')[0].strip('()')
        num, den = map(int, frac_str.split('/'))
        
        energy_correct = (num == expected_energy_numerator and den == expected_energy_denominator)
    except (ValueError, IndexError):
        # This handles cases where the energy string format is unexpected,
        # like in options C and D which have a different formula.
        return f"Incorrect energy formula. The energy for a 3D isotropic harmonic oscillator is proportional to (n + 3/2)hbar*omega, not '{answer_energy_str}'."

    # Verify the degeneracy value
    degeneracy_correct = (answer_degeneracy == expected_degeneracy)

    # --- Final Verdict ---
    if energy_correct and degeneracy_correct:
        return "Correct"
    else:
        error_messages = []
        if not energy_correct:
            error_messages.append(
                f"The energy value is incorrect. For the third excited state (n=3), the expected energy is ({expected_energy_numerator}/{expected_energy_denominator}) hbar*omega, but the answer gives a factor of ({num}/{den})."
            )
        if not degeneracy_correct:
            error_messages.append(
                f"The degeneracy is incorrect. For the n=3 state, the expected degeneracy is {expected_degeneracy}, but the answer gives {answer_degeneracy}."
            )
        return "\n".join(error_messages)

# Run the check and print the result
result = check_3d_harmonic_oscillator_answer()
print(result)