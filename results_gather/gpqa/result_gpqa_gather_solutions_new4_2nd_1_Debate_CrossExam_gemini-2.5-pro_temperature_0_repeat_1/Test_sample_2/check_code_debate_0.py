import numpy as np

def check_quantum_expectation_value():
    """
    This function checks the correctness of the calculated expectation value for a given quantum mechanics problem.

    The problem is to find the expectation value of the operator O = 10*sigma_z + 5*sigma_x
    for the state |ψ⟩ = 0.5|↑⟩ + sqrt(3)/2|↓⟩.

    The function performs the calculation using matrix representations and compares the result
    to the provided answer.
    """
    # Define the state vector |ψ⟩ in the standard basis where |↑⟩ = [1, 0] and |↓⟩ = [0, 1]
    # |ψ⟩ = 0.5 * [1, 0] + (sqrt(3)/2) * [0, 1] = [0.5, sqrt(3)/2]
    c_up = 0.5
    c_down = np.sqrt(3) / 2
    psi = np.array([c_up, c_down])

    # Check if the state is normalized (|c_up|^2 + |c_down|^2 = 1)
    norm = np.linalg.norm(psi)
    if not np.isclose(norm, 1.0):
        return f"Constraint not satisfied: The state vector is not normalized. Its norm is {norm:.4f}, but it should be 1."

    # Define the Pauli matrices
    sigma_z = np.array([[1, 0], [0, -1]])
    sigma_x = np.array([[0, 1], [1, 0]])

    # Construct the full operator O = 10*sigma_z + 5*sigma_x
    operator_O = 10 * sigma_z + 5 * sigma_x

    # Calculate the expectation value ⟨ψ|O|ψ⟩
    # This is calculated as psi_dagger * O * psi. Since psi is real, psi_dagger is just psi.T
    expectation_value = psi.conj().T @ operator_O @ psi
    
    # The question asks for the value up to one decimal place
    rounded_expectation_value = round(expectation_value, 1)

    # The provided answer is 'C', which corresponds to the value -0.7 from the options.
    # Options: A) 0.85, B) 1.65, C) -0.7, D) -1.4
    correct_answer_value = -0.7

    # Check if the calculated value matches the answer's value
    if np.isclose(rounded_expectation_value, correct_answer_value):
        return "Correct"
    else:
        # Provide a detailed reason for the discrepancy
        # Alternative calculation for clarity
        exp_val_z = c_up**2 - c_down**2
        exp_val_x = 2 * c_up * c_down
        calc_val = 10 * exp_val_z + 5 * exp_val_x
        
        reason = (f"Incorrect. The provided answer is {correct_answer_value}, but the calculated value is different.\n"
                  f"Step-by-step calculation:\n"
                  f"1. Expectation value of σz: ⟨σz⟩ = |0.5|² - |√3/2|² = 0.25 - 0.75 = -0.5\n"
                  f"2. Expectation value of σx: ⟨σx⟩ = ⟨ψ|σx|ψ⟩ = 2 * Re(c_up* * c_down) = 2 * 0.5 * (√3/2) = √3/2 ≈ {np.sqrt(3)/2:.4f}\n"
                  f"3. Total expectation value: ⟨O⟩ = 10*⟨σz⟩ + 5*⟨σx⟩ = 10*(-0.5) + 5*(√3/2) = -5 + 2.5*√3\n"
                  f"4. Numerical value: -5 + 2.5 * {np.sqrt(3):.4f} ≈ {calc_val:.4f}\n"
                  f"5. Rounding to one decimal place gives {rounded_expectation_value}.\n"
                  f"This does not match the answer's value of {correct_answer_value}.")
        return reason

# Run the check
result = check_quantum_expectation_value()
print(result)