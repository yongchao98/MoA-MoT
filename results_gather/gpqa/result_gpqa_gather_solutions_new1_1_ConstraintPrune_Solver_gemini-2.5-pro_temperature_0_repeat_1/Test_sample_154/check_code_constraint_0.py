import numpy as np

def check_correctness():
    """
    This function checks the correctness of the provided answer by recalculating the uncertainty of P_z.
    """
    # Let hbar = 1 for simplicity in calculation. The final answer will be a multiple of hbar.
    hbar = 1.0

    # --- Step 1: Define the problem's constants from the question ---
    # The state vector |psi> is given.
    psi = np.array([[-1/2], [1/np.sqrt(2)], [-1/2]])

    # The operator P_z is given.
    Pz = hbar * np.array([[1, 0, 0],
                          [0, 0, 0],
                          [0, 0, -1]])

    # The multiple-choice options from the question.
    # Note: The question text in the prompt has options A, B, C, D.
    # We will use these as the ground truth for mapping.
    options = {
        'A': hbar / np.sqrt(2),
        'B': np.sqrt(2) * hbar,
        'C': hbar / 2,
        'D': hbar
    }
    
    # The final answer provided by the LLM to be checked is 'A'.
    llm_answer_key = 'A'

    # --- Step 2: Perform the calculation based on quantum mechanics principles ---

    # Constraint: The state vector must be normalized, i.e., <psi|psi> = 1.
    norm_squared = (psi.conj().T @ psi)[0, 0]
    if not np.isclose(norm_squared, 1.0):
        return f"Incorrect: The state vector |psi> is not normalized. The inner product <psi|psi> is {norm_squared}, but it should be 1."

    # Constraint: Calculate the expectation value of Pz, <Pz> = <psi|Pz|psi>.
    exp_Pz = (psi.conj().T @ Pz @ psi)[0, 0]
    
    # The analysis correctly states <Pz> should be 0.
    if not np.isclose(exp_Pz, 0.0):
        return f"Incorrect: The calculation of the expectation value <Pz> is wrong. The calculated value is {exp_Pz}, but it should be 0."

    # Constraint: Calculate the expectation value of Pz^2, <Pz^2> = <psi|Pz^2|psi>.
    Pz_squared = Pz @ Pz
    exp_Pz2 = (psi.conj().T @ Pz_squared @ psi)[0, 0]

    # The analysis correctly states <Pz^2> should be hbar^2 / 2.
    if not np.isclose(exp_Pz2, hbar**2 / 2):
        return f"Incorrect: The calculation of the expectation value <Pz^2> is wrong. The calculated value is {exp_Pz2}, but it should be {hbar**2 / 2}."

    # Constraint: Calculate the uncertainty Delta Pz = sqrt(<Pz^2> - <Pz>^2).
    variance = exp_Pz2 - exp_Pz**2
    # Handle potential floating point inaccuracies resulting in a small negative number before sqrt
    if variance < 0 and np.isclose(variance, 0):
        variance = 0
    
    if variance < 0:
        return f"Incorrect: Calculation resulted in a negative variance ({variance}), which is physically impossible."
        
    uncertainty = np.sqrt(variance)

    # The analysis correctly states Delta Pz should be hbar / sqrt(2).
    if not np.isclose(uncertainty, hbar / np.sqrt(2)):
        return f"Incorrect: The final uncertainty calculation is wrong. The calculated value is {uncertainty}, but it should be {hbar / np.sqrt(2)}."

    # --- Step 3: Verify the final answer against the options ---
    
    # Find which option key corresponds to the calculated uncertainty.
    correct_option_key = None
    for key, value in options.items():
        if np.isclose(uncertainty, value):
            correct_option_key = key
            break
    
    if correct_option_key is None:
        return f"Incorrect: The calculated uncertainty {uncertainty}*hbar does not match any of the provided options."

    # Check if the LLM's chosen option matches the correct option.
    if llm_answer_key == correct_option_key:
        return "Correct"
    else:
        return f"Incorrect: The final answer choice is wrong. The calculated uncertainty is {uncertainty:.4f}*hbar, which corresponds to option {correct_option_key}. The provided answer was {llm_answer_key}."

# Run the check
result = check_correctness()
print(result)