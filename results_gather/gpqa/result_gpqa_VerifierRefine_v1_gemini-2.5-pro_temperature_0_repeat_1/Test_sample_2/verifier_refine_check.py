import numpy as np
import math

def check_correctness():
    """
    This function checks the correctness of the LLM's answer to the quantum mechanics problem.
    
    The problem is to find the expectation value of the operator 10*sigma_z + 5*sigma_x
    for a spin-half particle in the state |psi> = 0.5|up> + sqrt(3)/2|down>.
    """
    
    # --- Step 1: Define the quantum state and check normalization ---
    # The state is |psi> = 0.5|up> + sqrt(3)/2|down>
    # In the standard basis, |up> = [1, 0]^T and |down> = [0, 1]^T
    c1 = 0.5
    c2 = math.sqrt(3) / 2
    
    # The state vector |psi>
    psi = np.array([[c1], [c2]], dtype=complex)
    
    # The bra vector <psi| is the conjugate transpose of |psi>
    psi_bra = psi.conj().T
    
    # Constraint Check: The state must be normalized, i.e., <psi|psi> = 1
    # This means the sum of the squares of the absolute values of the coefficients must be 1.
    norm = np.vdot(psi, psi).real
    if not np.isclose(norm, 1.0):
        return f"Constraint not satisfied: The state vector is not normalized. The sum of the squares of the coefficients is {norm}, but it should be 1."

    # --- Step 2: Define the operator ---
    # The operator is O = 10*sigma_z + 5*sigma_x
    # Define the Pauli matrices
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    
    # Construct the full operator O
    operator_O = 10 * sigma_z + 5 * sigma_x

    # --- Step 3: Calculate the expectation value ---
    # The expectation value is <O> = <psi|O|psi>
    expectation_value_matrix = psi_bra @ operator_O @ psi
    
    # The result is a 1x1 matrix, so we extract the scalar value.
    # The expectation value of a Hermitian operator is always real.
    calculated_value = expectation_value_matrix[0, 0].real

    # --- Step 4: Compare with the LLM's answer ---
    # The LLM's answer corresponds to option B, which is -0.7.
    llm_answer = -0.7
    
    # Constraint Check: The question asks for the answer up to one decimal place.
    # We need to round our calculated value to one decimal place.
    rounded_calculated_value = round(calculated_value, 1)
    
    # Check if the rounded calculated value matches the LLM's answer.
    if np.isclose(rounded_calculated_value, llm_answer):
        return "Correct"
    else:
        exact_value_expr = f"-5 + 5*sqrt(3)/2"
        return (f"Incorrect. The calculated expectation value is {exact_value_expr} ≈ {calculated_value:.4f}. "
                f"Rounding this to one decimal place gives {rounded_calculated_value}. "
                f"The provided answer from option B is {llm_answer}, which does not match the rounded calculated value.")

# Execute the check and print the result
result = check_correctness()
print(result)