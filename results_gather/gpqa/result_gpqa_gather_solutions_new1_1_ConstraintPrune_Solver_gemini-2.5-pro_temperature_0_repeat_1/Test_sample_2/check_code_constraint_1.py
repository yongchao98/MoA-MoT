import numpy as np

def check_correctness():
    """
    This function checks the correctness of the answer to the quantum mechanics problem.

    The problem asks for the expectation value of the operator 10*sigma_z + 5*sigma_x
    for a spin-half particle in the state |ψ⟩ = 0.5|↑⟩ + sqrt(3)/2|↓⟩.
    The result should be rounded to one decimal place.
    """
    
    # Define the state vector |ψ⟩ in the z-basis where |↑⟩ = [1, 0] and |↓⟩ = [0, 1]
    # |ψ⟩ = 0.5 * [1, 0] + (sqrt(3)/2) * [0, 1]
    psi = np.array([0.5, np.sqrt(3)/2], dtype=complex)

    # Define the Pauli matrices
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)

    # Construct the full operator O = 10*sigma_z + 5*sigma_x
    operator_O = 10 * sigma_z + 5 * sigma_x

    # The expectation value ⟨O⟩ is calculated as ⟨ψ|O|ψ⟩.
    # ⟨ψ| is the conjugate transpose of |ψ⟩.
    bra_psi = psi.conj().T
    
    # Perform the matrix multiplication: ⟨ψ| * O * |ψ⟩
    expectation_value = bra_psi @ operator_O @ psi
    
    # The result should be a real number
    if not np.isclose(expectation_value.imag, 0):
        return f"Calculation error: The expectation value of a Hermitian operator must be real, but got {expectation_value}."

    real_expectation_value = expectation_value.real

    # Round the result to one decimal place as requested by the question
    rounded_value = round(real_expectation_value, 1)

    # The options provided in the question
    options = {
        'A': 0.85,
        'B': -1.4,
        'C': -0.7,
        'D': 1.65
    }

    # The final answer provided by the LLM is 'C'
    llm_answer_label = 'C'
    llm_answer_value = options[llm_answer_label]

    # Check if the calculated rounded value matches the value of the chosen option
    if np.isclose(rounded_value, llm_answer_value):
        return "Correct"
    else:
        # Find which option the calculation actually corresponds to
        correct_label = None
        for label, value in options.items():
            if np.isclose(rounded_value, value):
                correct_label = label
                break
        
        if correct_label:
            return (f"Incorrect. The calculated expectation value is {real_expectation_value:.4f}, "
                    f"which rounds to {rounded_value}. This corresponds to option {correct_label}. "
                    f"The provided answer was {llm_answer_label}.")
        else:
            return (f"Incorrect. The calculated expectation value is {real_expectation_value:.4f}, "
                    f"which rounds to {rounded_value}. This value does not match any of the options A, B, C, or D. "
                    f"The provided answer was {llm_answer_label}.")

# Run the check and print the result
result = check_correctness()
print(result)