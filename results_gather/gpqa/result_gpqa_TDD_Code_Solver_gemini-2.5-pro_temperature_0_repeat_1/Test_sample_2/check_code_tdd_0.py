import numpy as np

def check_quantum_expectation_value():
    """
    Checks the correctness of the calculated expectation value for the given quantum state and operator.
    """
    # The answer from the LLM is B, which corresponds to the value -0.7.
    llm_answer_value = -0.7

    # --- Step 1: Define the quantum state vector ---
    # The state is |ψ⟩ = 0.5|↑⟩ + sqrt(3)/2|↓⟩
    # In the standard basis where |↑⟩ = [1, 0] and |↓⟩ = [0, 1],
    # the state vector is [0.5, sqrt(3)/2].
    psi = np.array([0.5, np.sqrt(3)/2])

    # --- Step 2: Check if the state is normalized ---
    # The sum of the squares of the amplitudes must be 1.
    # (0.5)^2 + (sqrt(3)/2)^2 = 0.25 + 3/4 = 1.0
    norm = np.linalg.norm(psi)
    if not np.isclose(norm, 1.0):
        return f"Incorrect. The state vector |ψ⟩ is not normalized. Its norm is {norm:.4f}, but it should be 1."

    # --- Step 3: Define the operators as matrices ---
    sigma_z = np.array([[1, 0], [0, -1]])
    sigma_x = np.array([[0, 1], [1, 0]])

    # The full operator is O = 10σ_z + 5σ_x
    operator_O = 10 * sigma_z + 5 * sigma_x

    # --- Step 4: Calculate the expectation value ⟨ψ|O|ψ⟩ ---
    # ⟨ψ| is the conjugate transpose of |ψ⟩. Since the components are real, it's just the transpose.
    psi_bra = psi.conj().T
    
    # The expectation value is the matrix product: ψ_bra * O * ψ
    expectation_value = psi_bra @ operator_O @ psi
    
    # The question asks for the value up to one decimal place.
    calculated_value_rounded = round(expectation_value, 1)

    # --- Step 5: Compare the result with the LLM's answer ---
    if np.isclose(calculated_value_rounded, llm_answer_value):
        return "Correct"
    else:
        # Provide a detailed reason for the error.
        # Let's calculate the exact value symbolically for clarity:
        # ⟨O⟩ = 10 * ( |0.5|^2 - |sqrt(3)/2|^2 ) + 5 * (2 * Re(0.5 * sqrt(3)/2))
        # ⟨O⟩ = 10 * (0.25 - 0.75) + 5 * (sqrt(3)/2)
        # ⟨O⟩ = 10 * (-0.5) + 2.5 * sqrt(3)
        # ⟨O⟩ = -5 + 2.5 * sqrt(3)
        exact_value = -5 + 2.5 * np.sqrt(3)
        
        return (f"Incorrect. The calculated expectation value is {expectation_value:.4f} "
                f"(exact value: -5 + 2.5√3), which rounds to {calculated_value_rounded:.1f}. "
                f"The provided answer option B corresponds to {llm_answer_value}, which does not match the rounded calculation.")

# Run the check
result = check_quantum_expectation_value()
print(result)