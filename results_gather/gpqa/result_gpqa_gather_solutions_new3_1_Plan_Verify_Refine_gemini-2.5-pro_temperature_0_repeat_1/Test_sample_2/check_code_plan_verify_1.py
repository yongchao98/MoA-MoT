import numpy as np

def check_quantum_expectation_value():
    """
    Calculates the expectation value for the given quantum mechanics problem
    and checks if the provided answer is correct.
    """
    # --- Problem Definition ---
    # State |ψ⟩ = 0.5|↑⟩ + (sqrt(3)/2)|↓⟩
    c_up = 0.5
    c_down = np.sqrt(3) / 2
    
    # In the z-basis, |↑⟩ = [1, 0] and |↓⟩ = [0, 1]
    psi_ket = np.array([c_up, c_down], dtype=complex)

    # Pauli matrices
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)

    # The operator O = 10*σ_z + 5*σ_x
    operator_O = 10 * sigma_z + 5 * sigma_x

    # The options provided in the question
    options = {
        "A": 0.85,
        "B": -0.7,
        "C": -1.4,
        "D": 1.65
    }
    
    # The final answer from the LLM to be checked
    llm_answer_letter = "B"

    # --- Calculation ---
    # First, check if the state is normalized
    norm = np.linalg.norm(psi_ket)
    if not np.isclose(norm, 1.0):
        return f"Constraint failed: The state vector is not normalized. Norm is {norm:.4f}."

    # Bra vector ⟨ψ| is the conjugate transpose of the ket |ψ⟩
    psi_bra = psi_ket.conj().T
    
    # Expectation value is ⟨ψ|O|ψ⟩. The result must be real for a Hermitian operator.
    expectation_value = (psi_bra @ operator_O @ psi_ket).real

    # The question requires rounding to one decimal place
    rounded_value = round(expectation_value, 1)

    # --- Verification ---
    # Find which option letter corresponds to our calculated value
    correct_option_letter = None
    for letter, value in options.items():
        if np.isclose(rounded_value, value):
            correct_option_letter = letter
            break
    
    if correct_option_letter is None:
        return f"Calculation Error: The calculated value {rounded_value} does not match any of the provided options."

    # Check if the LLM's answer matches the correct option
    if llm_answer_letter == correct_option_letter:
        return "Correct"
    else:
        return (f"Incorrect. The calculated expectation value is {expectation_value:.4f}, "
                f"which rounds to {rounded_value}. This corresponds to option {correct_option_letter} "
                f"({options[correct_option_letter]}). The provided answer was {llm_answer_letter}.")

# Execute the check
result = check_quantum_expectation_value()
print(result)