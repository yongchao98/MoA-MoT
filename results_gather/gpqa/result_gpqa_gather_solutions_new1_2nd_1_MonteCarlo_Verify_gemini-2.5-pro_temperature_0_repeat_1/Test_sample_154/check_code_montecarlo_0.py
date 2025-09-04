import numpy as np

def check_correctness():
    """
    This function verifies the calculation of the uncertainty Delta Pz for the given quantum state.
    It follows these steps:
    1.  Defines the state vector |psi> and the operator Pz from the problem description.
    2.  Verifies that the state vector is normalized, a fundamental requirement for a physical state.
    3.  Verifies the consistency of the problem statement, i.e., that the given state is indeed an eigenstate of Px with the specified eigenvalue. This is not strictly necessary for the calculation but ensures the problem is well-posed.
    4.  Calculates the expectation value of Pz, <Pz>.
    5.  Calculates the expectation value of Pz^2, <Pz^2>.
    6.  Computes the uncertainty Delta Pz using the formula sqrt(<Pz^2> - <Pz>^2).
    7.  Compares the calculated numerical result with the value corresponding to the provided answer 'B'.
    8.  Returns "Correct" if the calculation matches the answer, or an error message otherwise.
    """
    # Set hbar = 1 for simplicity in numerical calculations, as all answers are in terms of hbar.
    hbar = 1.0

    # --- Define quantities from the problem statement ---

    # The state of the system, |psi>
    psi_ket = np.array([[-1/2.0], 
                        [1/np.sqrt(2)], 
                        [-1/2.0]])

    # The operator Pz
    Pz = hbar * np.array([[1.0, 0.0, 0.0],
                          [0.0, 0.0, 0.0],
                          [0.0, 0.0, -1.0]])

    # The options as listed in the final provided answer's analysis
    options = {
        'A': hbar / 2.0,
        'B': hbar / np.sqrt(2),
        'C': hbar,
        'D': np.sqrt(2) * hbar
    }
    
    # The final answer given by the LLM to be checked
    llm_answer_label = 'B'

    # --- Verification Steps ---

    # 1. Check for state normalization: <psi|psi> = 1
    psi_bra = psi_ket.conj().T
    norm_squared = (psi_bra @ psi_ket).item()
    if not np.isclose(norm_squared, 1.0):
        return f"Constraint not satisfied: The state vector |psi> is not normalized. <psi|psi> = {norm_squared:.4f}, but it should be 1."

    # 2. Check for problem consistency (optional but good practice)
    # The problem states |psi> is an eigenstate of Px with eigenvalue -hbar.
    Px = (hbar / np.sqrt(2)) * np.array([[0.0, 1.0, 0.0],
                                         [1.0, 0.0, 1.0],
                                         [0.0, 1.0, 0.0]])
    # Calculate Px|psi>
    px_psi = Px @ psi_ket
    # Calculate the expected result: -hbar * |psi>
    expected_eigenstate_result = -hbar * psi_ket
    if not np.allclose(px_psi, expected_eigenstate_result):
        return f"Constraint not satisfied: The problem statement is inconsistent. The given state |psi> is not an eigenstate of Px with eigenvalue -hbar as claimed."

    # 3. Calculate the expectation value <Pz>
    exp_Pz = (psi_bra @ Pz @ psi_ket).item()

    # 4. Calculate the expectation value <Pz^2>
    Pz_squared = Pz @ Pz
    exp_Pz_squared = (psi_bra @ Pz_squared @ psi_ket).item()

    # 5. Calculate the uncertainty Delta Pz
    # (Delta Pz)^2 = <Pz^2> - <Pz>^2
    variance_Pz = exp_Pz_squared - exp_Pz**2
    
    # Handle potential floating point issues where variance is a tiny negative number
    if variance_Pz < 0 and np.isclose(variance_Pz, 0):
        variance_Pz = 0
    
    if variance_Pz < 0:
        return f"Incorrect: Calculation resulted in a negative variance ({variance_Pz:.4f}), which is physically impossible."

    uncertainty_Pz = np.sqrt(variance_Pz)

    # 6. Compare the calculated uncertainty with the value from the chosen option
    expected_value = options.get(llm_answer_label)
    
    if expected_value is None:
        return f"Incorrect: The provided answer label '{llm_answer_label}' is not a valid option (A, B, C, D)."

    if np.isclose(uncertainty_Pz, expected_value):
        return "Correct"
    else:
        # Find which option the calculation actually matches
        correct_label = "None"
        for label, val in options.items():
            if np.isclose(uncertainty_Pz, val):
                correct_label = label
                break
        
        return (f"Incorrect: The calculated uncertainty is Delta_Pz = {uncertainty_Pz:.4f} * hbar. "
                f"This corresponds to option {correct_label} ({options.get(correct_label, 0):.4f} * hbar). "
                f"The provided answer was {llm_answer_label} ({expected_value:.4f} * hbar).")

# Run the check
result = check_correctness()
print(result)