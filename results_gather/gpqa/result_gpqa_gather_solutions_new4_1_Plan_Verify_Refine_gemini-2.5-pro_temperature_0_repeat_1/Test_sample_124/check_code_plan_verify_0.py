import sympy

def check_answer():
    """
    Checks the correctness of the answer for the quantum harmonic oscillator problem.
    """
    # Define the constants and variables from the problem
    # We don't need their numerical values, just their symbols for verification
    hbar = sympy.Symbol('hbar')
    omega = sympy.Symbol('omega')

    # The question asks for the third excited state.
    # Ground state: n = 0
    # First excited state: n = 1
    # Second excited state: n = 2
    # Third excited state: n = 3
    n = 3

    # --- Step 1: Calculate the correct energy ---
    # The energy eigenvalues for a 3D isotropic harmonic oscillator are given by:
    # E_n = (n + 3/2) * hbar * omega
    # For the third excited state (n=3):
    # E_3 = (3 + 3/2) * hbar * omega = (9/2) * hbar * omega
    correct_energy_expr = (sympy.Rational(9, 2)) * hbar * omega
    
    # --- Step 2: Calculate the correct degeneracy ---
    # The degeneracy for the n-th level of a 3D isotropic harmonic oscillator is:
    # g_n = (n + 1) * (n + 2) / 2
    # For the third excited state (n=3):
    # g_3 = (3 + 1) * (3 + 2) / 2 = (4 * 5) / 2 = 10
    correct_degeneracy = int(((n + 1) * (n + 2)) / 2)

    # --- Step 3: Define the options provided in the question ---
    # Note: The energy expressions are stored as strings for initial comparison
    options = {
        'A': ('11 \\pi^2 \\hbar^2 / (2m r^2)', 10),
        'B': ('11 \\pi^2 \\hbar^2 / (2m r^2)', 3),
        'C': ('(9/2) \\hbar \\omega', 10),
        'D': ('(9/2) \\hbar \\omega', 3)
    }

    # The provided answer to check is 'C'
    llm_answer_key = 'C'
    
    # --- Step 4: Retrieve the values from the chosen answer ---
    if llm_answer_key not in options:
        return f"The answer key '{llm_answer_key}' is not a valid option (A, B, C, D)."
        
    llm_energy_str, llm_degeneracy = options[llm_answer_key]

    # --- Step 5: Verify the degeneracy ---
    if llm_degeneracy != correct_degeneracy:
        return (f"Incorrect. The degeneracy is wrong. "
                f"The calculated degeneracy for the third excited state (n=3) is {correct_degeneracy}, "
                f"but the answer provides a degeneracy of {llm_degeneracy}.")

    # --- Step 6: Verify the energy expression ---
    # Check if the energy formula corresponds to a harmonic oscillator.
    # The correct formula depends on hbar and omega, not 1/r^2.
    if "omega" not in llm_energy_str or "r^2" in llm_energy_str:
        return (f"Incorrect. The energy expression '{llm_energy_str}' is for the wrong physical system. "
                f"The potential is for a harmonic oscillator, so the energy should be proportional to hbar*omega, "
                f"not 1/r^2 (which is for a particle in a spherical well).")

    # Check the numerical coefficient in the energy expression.
    # We expect the coefficient to be 9/2 = 4.5.
    # A simple way to check is to see if '9/2' is in the string.
    if '9/2' not in llm_energy_str:
        return (f"Incorrect. The numerical coefficient in the energy expression is wrong. "
                f"Expected a coefficient of 9/2 for the third excited state, but the answer's expression is '{llm_energy_str}'.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)