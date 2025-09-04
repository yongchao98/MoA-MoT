import numpy as np

def check_quantum_expectation_value():
    """
    This function checks the correctness of the calculated expectation value for the given quantum mechanics problem.

    The problem is to find the expectation value of the operator O = 10σ_z + 5σ_x
    for the state |ψ⟩ = 0.5|↑⟩ + sqrt(3)/2|↓⟩.

    The final answer provided is 'B', which corresponds to the value -0.7.
    """

    # Define the options as provided in the question
    options = {
        'A': 0.85,
        'B': -0.7,
        'C': -1.4,
        'D': 1.65
    }

    # The final answer to be checked is 'B'
    provided_answer_key = 'B'
    provided_answer_value = options[provided_answer_key]

    try:
        # --- Step 1: Define the quantum state and operators in matrix form ---

        # Coefficients of the state vector
        c_up = 0.5
        c_down = np.sqrt(3) / 2

        # State vector |ψ⟩ in the z-basis where |↑⟩ = [1, 0]^T and |↓⟩ = [0, 1]^T
        psi_ket = np.array([[c_up], [c_down]], dtype=complex)

        # Bra vector ⟨ψ| (conjugate transpose of |ψ⟩)
        psi_bra = psi_ket.conj().T

        # Check for normalization (the sum of squared magnitudes of coefficients must be 1)
        norm_squared = np.vdot(psi_ket, psi_ket).real
        if not np.isclose(norm_squared, 1.0):
            return f"Constraint not satisfied: The state vector is not normalized. The sum of the squared coefficients is {norm_squared:.4f}, but it should be 1."

        # Pauli matrices
        sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
        sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)

        # --- Step 2: Calculate the expectation value ---

        # Expectation value of σ_z: ⟨σ_z⟩ = ⟨ψ|σ_z|ψ⟩
        exp_val_z = (psi_bra @ sigma_z @ psi_ket)[0, 0].real

        # Expectation value of σ_x: ⟨σ_x⟩ = ⟨ψ|σ_x|ψ⟩
        exp_val_x = (psi_bra @ sigma_x @ psi_ket)[0, 0].real

        # Total expectation value of O = 10σ_z + 5σ_x
        # Using linearity: ⟨O⟩ = 10⟨σ_z⟩ + 5⟨σ_x⟩
        total_expectation_value = 10 * exp_val_z + 5 * exp_val_x

        # --- Step 3: Compare the result with the provided answer ---

        # The question asks for the value up to one decimal place
        calculated_rounded_value = round(total_expectation_value, 1)

        # Check if the calculated rounded value matches the value from the chosen option 'B'
        if np.isclose(calculated_rounded_value, provided_answer_value):
            return "Correct"
        else:
            return (f"Incorrect. The calculated expectation value is {total_expectation_value:.4f}, "
                    f"which rounds to {calculated_rounded_value}. The provided answer is '{provided_answer_key}', "
                    f"which corresponds to the value {provided_answer_value}. The calculated value does not match the answer's value.")

    except Exception as e:
        return f"An error occurred during the calculation: {e}"

# Run the check
result = check_quantum_expectation_value()
print(result)