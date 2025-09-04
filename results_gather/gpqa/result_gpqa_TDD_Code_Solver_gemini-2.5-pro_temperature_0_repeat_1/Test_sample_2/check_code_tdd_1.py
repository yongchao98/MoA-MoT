import numpy as np

def check_quantum_expectation_value():
    """
    Checks the correctness of the calculated expectation value for the given quantum state and operator.
    """
    # The state is |ψ⟩ = 0.5|↑⟩ + (sqrt(3)/2)|↓⟩
    # In the standard basis where |↑⟩ = [1, 0]^T and |↓⟩ = [0, 1]^T,
    # the state vector |ψ⟩ is:
    psi = np.array([[0.5], [np.sqrt(3)/2]])

    # First, let's check if the state is normalized, which is a fundamental constraint.
    # The norm squared should be ⟨ψ|ψ⟩ = 1.
    norm_squared = np.vdot(psi, psi).real
    if not np.isclose(norm_squared, 1.0):
        return f"Constraint failed: The given state |ψ⟩ is not normalized. Its norm squared is {norm_squared:.4f}, but it should be 1."

    # Define the Pauli matrices
    sigma_z = np.array([[1, 0], [0, -1]])
    sigma_x = np.array([[0, 1], [1, 0]])

    # Define the operator from the question: O = 10σ_z + 5σ_x
    operator = 10 * sigma_z + 5 * sigma_x

    # Calculate the expectation value ⟨ψ|O|ψ⟩
    # ⟨ψ| is the conjugate transpose of |ψ⟩
    psi_bra = psi.conj().T
    
    # The expectation value is calculated by the matrix product ⟨ψ| * O * |ψ⟩
    expectation_value_matrix = psi_bra @ operator @ psi
    
    # The result is a 1x1 matrix, so we extract the scalar value
    calculated_value = expectation_value_matrix[0, 0]

    # The question asks for the value up to one decimal place
    rounded_value = round(calculated_value, 1)

    # The provided answer is -0.7, which corresponds to option B.
    expected_value = -0.7

    # Check if the calculated value matches the answer's value.
    # We use np.isclose for robust floating-point comparison.
    if np.isclose(rounded_value, expected_value):
        # Let's also verify the intermediate steps in the provided answer's logic.
        # ⟨ψ|10σ_z|ψ⟩ = -5
        exp_10_sigma_z = (psi_bra @ (10 * sigma_z) @ psi)[0, 0]
        if not np.isclose(exp_10_sigma_z, -5.0):
            return f"Incorrect. The final answer is correct, but the intermediate step for ⟨ψ|10σ_z|ψ⟩ is wrong. Calculated: {exp_10_sigma_z}, Expected: -5.0"
        
        # ⟨ψ|5σ_x|ψ⟩ = 5 * sqrt(3)/2 ≈ 4.33
        exp_5_sigma_x = (psi_bra @ (5 * sigma_x) @ psi)[0, 0]
        if not np.isclose(exp_5_sigma_x, 5 * np.sqrt(3) / 2):
             return f"Incorrect. The final answer is correct, but the intermediate step for ⟨ψ|5σ_x|ψ⟩ is wrong. Calculated: {exp_5_sigma_x}, Expected: {5 * np.sqrt(3) / 2}"

        # Total expectation value = -5 + 5*sqrt(3)/2 ≈ -0.67
        total_unrounded = exp_10_sigma_z + exp_5_sigma_x
        if not np.isclose(calculated_value, total_unrounded):
            return f"Incorrect. The final answer is correct, but the unrounded total is inconsistent. Calculated: {calculated_value}, Expected from steps: {total_unrounded}"

        return "Correct"
    else:
        return f"Incorrect. The calculated expectation value is {calculated_value:.4f}, which rounds to {rounded_value}. This does not match the provided answer's value of {expected_value}."

# Execute the check
result = check_quantum_expectation_value()
print(result)