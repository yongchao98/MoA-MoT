import numpy as np

def check_quantum_expectation_value():
    """
    This function calculates the expectation value of the operator 10σ_z + 5σ_x
    for the given spin state and checks if the provided answer is correct.
    """
    # --- 1. Define the state and operators ---
    
    # The state is |ψ⟩ = 0.5|↑⟩ + (sqrt(3)/2)|↓⟩
    # In the standard basis where |↑⟩ = [1, 0]^T and |↓⟩ = [0, 1]^T,
    # the state vector is:
    psi_ket = np.array([0.5, np.sqrt(3)/2])
    
    # The bra vector ⟨ψ| is the conjugate transpose of the ket.
    # Since the components are real, it's just the transpose.
    psi_bra = psi_ket.T

    # The Pauli matrices
    sigma_z = np.array([[1, 0], [0, -1]])
    sigma_x = np.array([[0, 1], [1, 0]])

    # --- 2. Verify Constraints ---

    # Constraint: The state must be normalized (||ψ|| = 1)
    norm = np.linalg.norm(psi_ket)
    if not np.isclose(norm, 1.0):
        return f"Incorrect: The state vector is not normalized. Its norm is {norm:.4f}, but it should be 1."

    # --- 3. Perform the Calculation ---

    # Due to linearity, ⟨10σ_z + 5σ_x⟩ = 10⟨σ_z⟩ + 5⟨σ_x⟩
    
    # Expectation value of σ_z: ⟨ψ|σ_z|ψ⟩
    exp_sigma_z = psi_bra @ sigma_z @ psi_ket
    
    # Expectation value of σ_x: ⟨ψ|σ_x|ψ⟩
    exp_sigma_x = psi_bra @ sigma_x @ psi_ket
    
    # Total expectation value
    total_expectation_value = 10 * exp_sigma_z + 5 * exp_sigma_x
    
    # --- 4. Check the Final Answer ---

    # The question asks for the value up to one decimal place.
    calculated_rounded_value = round(total_expectation_value, 1)
    
    # The provided answer is 'A', which corresponds to -0.7.
    # Let's check if our calculation matches this.
    expected_value = -0.7
    
    if np.isclose(calculated_rounded_value, expected_value):
        return "Correct"
    else:
        # Provide a detailed reason for the failure.
        exact_value = -5 + 2.5 * np.sqrt(3)
        return (f"Incorrect: The provided answer is wrong.\n"
                f"The calculated expectation value of σ_z is {exp_sigma_z:.4f}.\n"
                f"The calculated expectation value of σ_x is {exp_sigma_x:.4f}.\n"
                f"The total exact expectation value is 10*({exp_sigma_z:.4f}) + 5*({exp_sigma_x:.4f}) = {total_expectation_value:.4f}.\n"
                f"This value rounded to one decimal place is {calculated_rounded_value}.\n"
                f"The provided answer corresponds to {expected_value}, which does not match the calculation.")

# Execute the check and print the result
result = check_quantum_expectation_value()
print(result)