import numpy as np

def check_quantum_expectation_value():
    """
    This function calculates the expectation value of the operator 10σ_z + 5σ_x
    for the given quantum state and checks if it matches the provided answer.
    """
    
    # 1. Define the quantum state |ψ⟩ = 0.5|↑⟩ + (√3)/2|↓⟩
    # In the standard basis where |↑⟩ = [1, 0] and |↓⟩ = [0, 1],
    # the state vector is:
    psi = np.array([0.5, np.sqrt(3)/2], dtype=complex)

    # Sanity check: Ensure the state is normalized (||ψ||^2 = 1)
    norm_squared = np.vdot(psi, psi).real
    if not np.isclose(norm_squared, 1.0):
        return f"Constraint not satisfied: The state vector is not normalized. Its norm squared is {norm_squared:.4f}."

    # 2. Define the Pauli matrices σ_z and σ_x
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)

    # 3. Construct the full operator O = 10σ_z + 5σ_x
    O = 10 * sigma_z + 5 * sigma_x

    # 4. Calculate the expectation value ⟨O⟩ = ⟨ψ|O|ψ⟩
    # ⟨ψ| is the conjugate transpose of |ψ⟩
    psi_bra = psi.conj().T
    
    # The expectation value is calculated by the matrix product: ⟨ψ| * O * |ψ⟩
    expectation_value = (psi_bra @ O @ psi).real

    # 5. Round the result to one decimal place as per the question's constraint
    calculated_answer = round(expectation_value, 1)

    # 6. Define the options and the provided answer
    options = {'A': 0.85, 'B': -1.4, 'C': 1.65, 'D': -0.7}
    llm_answer_choice = 'D'
    
    expected_value_from_answer = options.get(llm_answer_choice)

    # 7. Check if the calculated answer matches the provided answer
    if np.isclose(calculated_answer, expected_value_from_answer):
        return "Correct"
    else:
        return (f"Incorrect. The calculated expectation value is {expectation_value:.4f}, "
                f"which rounds to {calculated_answer}. The provided answer choice '{llm_answer_choice}' "
                f"corresponds to the value {expected_value_from_answer}, which does not match.")

# Run the check
result = check_quantum_expectation_value()
print(result)