import numpy as np

def check_expectation_value():
    """
    Calculates the expectation value for the given quantum mechanics problem
    and checks if it matches the provided answer.
    """
    # 1. Define the quantum state
    # The state is |ψ⟩ = 0.5|↑⟩ + (sqrt(3)/2)|↓⟩.
    # In the standard basis, |↑⟩ = [1, 0]^T and |↓⟩ = [0, 1]^T.
    c_up = 0.5
    c_down = np.sqrt(3) / 2
    
    # The state vector (ket) is represented as a column vector.
    psi = np.array([[c_up], [c_down]], dtype=complex)

    # Sanity check: ensure the state is normalized (i.e., ⟨ψ|ψ⟩ = 1).
    norm_squared = (psi.conj().T @ psi).item()
    if not np.isclose(norm_squared, 1.0):
        return f"State is not normalized. Sum of squared coefficients is {norm_squared:.4f}."

    # 2. Define the operators
    # The Pauli matrices σ_z and σ_x.
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)

    # The full operator is O = 10σ_z + 5σ_x.
    operator_O = 10 * sigma_z + 5 * sigma_x

    # 3. Calculate the expectation value
    # The expectation value is given by ⟨O⟩ = ⟨ψ|O|ψ⟩.
    # ⟨ψ| (bra) is the conjugate transpose of |ψ⟩ (ket).
    psi_bra = psi.conj().T
    expectation_value_exact = (psi_bra @ operator_O @ psi).item()

    # The question asks for the result rounded to one decimal place.
    # The result is real, so we take the real part before rounding.
    calculated_result = round(np.real(expectation_value_exact), 1)

    # 4. Check against the provided answer
    # The provided answer is <<<A>>>, which corresponds to -0.7.
    answer_to_check = -0.7

    print(f"Exact calculated expectation value: {np.real(expectation_value_exact):.5f}")
    print(f"Calculated value rounded to one decimal place: {calculated_result}")
    print(f"The value to check (for option A) is: {answer_to_check}")

    # Compare the calculated result with the answer to be checked.
    if np.isclose(calculated_result, answer_to_check):
        return "Correct"
    else:
        return (f"Incorrect. The calculated value rounds to {calculated_result}, "
                f"but the provided answer 'A' corresponds to {answer_to_check}.")

# Run the check
result = check_expectation_value()
# The final verdict is printed below the code block as requested.