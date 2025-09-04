import numpy as np

def check_quantum_expectation_value():
    """
    This function checks the correctness of the calculated expectation value for a given quantum state and operator.
    
    The problem is:
    State: |ψ⟩ = 0.5|↑⟩ + sqrt(3)/2|↓⟩
    Operator: O = 10σ_z + 5σ_x
    Goal: Find ⟨ψ|O|ψ⟩ rounded to one decimal place.
    """
    try:
        # --- 1. Define the quantum state and operators in the z-basis ---
        
        # |↑⟩ is represented as the vector [1, 0]
        # |↓⟩ is represented as the vector [0, 1]
        
        # The state |ψ⟩ = 0.5|↑⟩ + sqrt(3)/2|↓⟩
        # This is a column vector.
        psi = np.array([[0.5], [np.sqrt(3)/2]])
        
        # The bra ⟨ψ| is the conjugate transpose of |ψ⟩.
        # Since the components are real, it's just the transpose.
        bra_psi = psi.T

        # Check if the state is normalized. The sum of squared amplitudes should be 1.
        norm = np.vdot(psi, psi).real
        if not np.isclose(norm, 1.0):
            return f"Constraint not satisfied: The state vector |ψ⟩ is not normalized. The norm-squared is {norm}, but it should be 1."

        # Define the Pauli matrices
        sigma_z = np.array([[1, 0], [0, -1]])
        sigma_x = np.array([[0, 1], [1, 0]])
        
        # --- 2. Construct the full operator ---
        
        # O = 10σ_z + 5σ_x
        operator_O = 10 * sigma_z + 5 * sigma_x
        
        # --- 3. Calculate the expectation value ⟨ψ|O|ψ⟩ ---
        
        # The expectation value is calculated as ⟨ψ|O|ψ⟩ = bra_psi * operator_O * psi
        expectation_value_matrix = bra_psi @ operator_O @ psi
        
        # The result is a 1x1 matrix, so we extract the scalar value.
        calculated_value = expectation_value_matrix[0, 0]
        
        # --- 4. Compare with the given answer ---
        
        # The given answer is A, which corresponds to the value -0.7.
        # The question asks for the value up to one decimal place.
        expected_value = -0.7
        
        # Round the calculated value to one decimal place for comparison.
        rounded_calculated_value = round(calculated_value, 1)
        
        # Check if the rounded calculated value matches the expected value.
        if np.isclose(rounded_calculated_value, expected_value):
            return "Correct"
        else:
            return (f"The answer is incorrect.\n"
                    f"The calculated expectation value is {calculated_value:.4f}.\n"
                    f"When rounded to one decimal place, this is {rounded_calculated_value}.\n"
                    f"The value from the selected answer 'A' is {expected_value}.\n"
                    f"The calculated value does not match the answer.")

    except Exception as e:
        return f"An error occurred during the calculation: {e}"

# Run the check
result = check_quantum_expectation_value()
print(result)