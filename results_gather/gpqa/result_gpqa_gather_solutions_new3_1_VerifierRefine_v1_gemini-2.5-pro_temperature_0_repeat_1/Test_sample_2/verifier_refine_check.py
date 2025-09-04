import numpy as np

def check_quantum_expectation_value():
    """
    This function calculates the expectation value for the given quantum mechanics problem
    and checks if the provided answer is correct.

    Problem:
    A spin-half particle is in a linear superposition |ψ⟩ = 0.5|↑⟩ + sqrt(3)/2|↓⟩.
    Calculate the expectation value of the operator O = 10σ_z + 5σ_x, rounded to one decimal place.

    Options:
    A) 0.85
    B) -0.7
    C) -1.4
    D) 1.65

    The provided final answer is <<<B>>>.
    """
    try:
        # 1. Define the state vector |ψ⟩
        # In the standard basis where |↑⟩ = [1, 0] and |↓⟩ = [0, 1]
        c_up = 0.5
        c_down = np.sqrt(3) / 2
        psi = np.array([c_up, c_down], dtype=complex)

        # 2. Check if the state is normalized (sum of squared magnitudes of coefficients is 1)
        norm = np.vdot(psi, psi).real
        if not np.isclose(norm, 1.0):
            return f"Constraint not satisfied: The state |ψ⟩ is not normalized. The sum of the squared coefficients is {norm:.4f}, but it should be 1."

        # 3. Define the Pauli matrices σ_z and σ_x
        sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
        sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)

        # 4. Calculate the expectation values for σ_z and σ_x separately
        # ⟨σ_z⟩ = ⟨ψ|σ_z|ψ⟩
        exp_val_z = np.vdot(psi, sigma_z @ psi).real
        
        # ⟨σ_x⟩ = ⟨ψ|σ_x|ψ⟩
        exp_val_x = np.vdot(psi, sigma_x @ psi).real

        # 5. Calculate the total expectation value for O = 10σ_z + 5σ_x
        # ⟨O⟩ = 10⟨σ_z⟩ + 5⟨σ_x⟩
        total_exp_val = 10 * exp_val_z + 5 * exp_val_x

        # 6. Round the result to one decimal place as required by the question
        rounded_result = round(total_exp_val, 1)

        # 7. Check against the provided answer
        # The provided answer is 'B', which corresponds to -0.7
        correct_value = -0.7
        
        if np.isclose(rounded_result, correct_value):
            return "Correct"
        else:
            # Provide a detailed reason for the error
            reason = (f"The final answer is incorrect.\n"
                      f"The calculated expectation value of σ_z is {exp_val_z:.4f}.\n"
                      f"The calculated expectation value of σ_x is {exp_val_x:.4f} (which is sqrt(3)/2).\n"
                      f"The total expectation value is 10*({exp_val_z:.4f}) + 5*({exp_val_x:.4f}) = {total_exp_val:.4f}.\n"
                      f"Rounding to one decimal place gives {rounded_result}.\n"
                      f"The provided answer corresponds to {correct_value}, which does not match the calculated value.")
            return reason

    except Exception as e:
        return f"An error occurred during the calculation: {e}"

# Execute the check and print the result
print(check_quantum_expectation_value())