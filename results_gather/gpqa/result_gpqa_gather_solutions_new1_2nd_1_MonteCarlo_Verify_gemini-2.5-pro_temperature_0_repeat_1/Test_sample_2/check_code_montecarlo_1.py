import numpy as np
import math

def check_correctness_of_quantum_expectation():
    """
    This function checks the correctness of the final answer for the given quantum mechanics problem.
    
    Question: A spin-half particle is in a linear superposition 0.5|u>+sqrt(3)/2|d> of its spin-up and spin-down states. 
    If |u> and |d> are the eigenstates of sigma_z, then what is the expectation value up to one decimal place, 
    of the operator 10*sigma_z + 5*sigma_x?
    
    Options: A) -1.4, B) -0.7, C) 0.85, D) 1.65
    Provided Answer: B
    """

    # --- Step 1: Define the state vector and Pauli matrices ---
    # The state is |ψ⟩ = 0.5|↑⟩ + (√3/2)|↓⟩
    # In the z-basis, |↑⟩ = [1, 0] and |↓⟩ = [0, 1]
    psi = np.array([0.5, math.sqrt(3)/2], dtype=complex)

    # Pauli matrices
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)

    # --- Step 2: Define the operator ---
    # The operator is O = 10*σ_z + 5*σ_x
    operator_O = 10 * sigma_z + 5 * sigma_x

    # --- Step 3: Check state normalization (Constraint) ---
    norm = np.vdot(psi, psi).real
    if not np.isclose(norm, 1.0):
        return f"Constraint not satisfied: The state vector is not normalized. Its norm is {norm:.4f}."

    # --- Step 4: Calculate the expectation value ⟨ψ|O|ψ⟩ ---
    # The bra vector ⟨ψ| is the conjugate transpose of the ket |ψ⟩
    psi_bra = np.conjugate(psi).T
    
    # Expectation value is ⟨ψ|O|ψ⟩
    expectation_value = (psi_bra @ operator_O @ psi).real

    # --- Step 5: Round to one decimal place (Constraint) ---
    calculated_rounded_value = round(expectation_value, 1)

    # --- Step 6: Compare with the provided answer ---
    # The provided answer is 'B', which corresponds to -0.7
    provided_answer_label = 'B'
    options = {'A': -1.4, 'B': -0.7, 'C': 0.85, 'D': 1.65}
    
    if provided_answer_label not in options:
        return f"Error: The provided answer label '{provided_answer_label}' is not a valid option."
        
    target_value = options[provided_answer_label]

    # Check if the calculated value matches the target value
    if np.isclose(calculated_rounded_value, target_value):
        return "Correct"
    else:
        reason = (
            f"The final answer is incorrect.\n"
            f"1. The calculated expectation value is {expectation_value:.5f}.\n"
            f"2. The question requires rounding to one decimal place, which gives {calculated_rounded_value}.\n"
            f"3. The provided answer is '{provided_answer_label}', which corresponds to the value {target_value}.\n"
            f"4. The calculated value {calculated_rounded_value} does not match the answer's value {target_value}."
        )
        return reason

# Execute the check
# print(check_correctness_of_quantum_expectation())