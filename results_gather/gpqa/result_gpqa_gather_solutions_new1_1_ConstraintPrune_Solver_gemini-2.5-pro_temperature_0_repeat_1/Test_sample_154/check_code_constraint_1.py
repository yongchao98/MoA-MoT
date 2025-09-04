import numpy as np

def check_answer():
    """
    Checks the correctness of the provided answer for the quantum mechanics problem.
    """
    # Let hbar = 1.0 for simplicity in numerical calculations.
    # The final result will be in units of hbar.
    hbar = 1.0

    # Constraint 1: Define the physical system from the question.
    # The state vector |psi>
    psi = np.array([[-1/2], [1/np.sqrt(2)], [-1/2]])
    # The operator P_z
    Pz = hbar * np.array([[1, 0, 0], [0, 0, 0], [0, 0, -1]])

    # The multiple-choice options given in the question analysis.
    options = {
        "A": hbar / np.sqrt(2),
        "B": np.sqrt(2) * hbar,
        "C": hbar / 2,
        "D": hbar
    }
    
    # The final answer provided by the LLM to be checked.
    llm_answer_key = "A"

    # --- Verification Steps ---

    # Step 1: Verify the state vector is normalized, a fundamental requirement for a physical state.
    # <psi|psi> should be 1.
    psi_bra = psi.conj().T
    norm_squared = (psi_bra @ psi)[0, 0]
    if not np.isclose(norm_squared, 1.0):
        return f"Incorrect. The state vector |psi> is not normalized. <psi|psi> = {norm_squared}, but it should be 1."

    # Step 2: Calculate the expectation value of Pz, <Pz>.
    # <Pz> = <psi|Pz|psi>
    exp_Pz = (psi_bra @ Pz @ psi)[0, 0]
    if not np.isclose(exp_Pz, 0.0):
        return f"Incorrect. The calculation of <Pz> is wrong. Expected 0, but calculated {exp_Pz}."

    # Step 3: Calculate the expectation value of Pz^2, <Pz^2>.
    # <Pz^2> = <psi|Pz^2|psi>
    Pz2 = Pz @ Pz
    exp_Pz2 = (psi_bra @ Pz2 @ psi)[0, 0]
    expected_exp_Pz2 = hbar**2 / 2
    if not np.isclose(exp_Pz2, expected_exp_Pz2):
        return f"Incorrect. The calculation of <Pz^2> is wrong. Expected {expected_exp_Pz2}*hbar^2, but calculated {exp_Pz2}*hbar^2."

    # Step 4: Calculate the uncertainty Delta Pz using the formula.
    # Delta Pz = sqrt(<Pz^2> - <Pz>^2)
    uncertainty_squared = exp_Pz2 - exp_Pz**2
    calculated_uncertainty = np.sqrt(uncertainty_squared)

    # Step 5: Compare the calculated uncertainty with the value corresponding to the LLM's answer key.
    llm_answer_value = options.get(llm_answer_key)
    if llm_answer_value is None:
        return f"Incorrect. The provided answer key '{llm_answer_key}' is not a valid option (A, B, C, D)."

    if np.isclose(calculated_uncertainty, llm_answer_value):
        # The calculation is correct and the mapping to the option is also correct.
        return "Correct"
    else:
        # Find which option the calculation actually matches.
        correct_key = "None"
        for key, value in options.items():
            if np.isclose(calculated_uncertainty, value):
                correct_key = key
                break
        
        return (f"Incorrect. The final answer is given as {llm_answer_key}, which corresponds to a value of {llm_answer_value:.4f}*hbar. "
                f"However, the correct calculated uncertainty is {calculated_uncertainty:.4f}*hbar, which corresponds to option {correct_key}.")

# Run the check
result = check_answer()
print(result)