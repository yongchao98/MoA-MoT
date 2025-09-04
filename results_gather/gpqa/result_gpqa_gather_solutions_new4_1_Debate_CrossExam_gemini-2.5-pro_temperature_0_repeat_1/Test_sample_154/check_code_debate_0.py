import numpy as np
import sympy

def check_correctness_of_answer():
    """
    Checks the correctness of the LLM's answer to the quantum mechanics problem.
    The code recalculates the uncertainty from the problem's premises and compares
    it to the provided answer.
    """
    # --- Problem Setup ---
    # Define hbar as a symbolic variable for exact mathematical representation
    hbar = sympy.Symbol('hbar')

    # Define the operator P_z matrix as given in the question
    Pz = sympy.Matrix([
        [hbar, 0, 0],
        [0, 0, 0],
        [0, 0, -hbar]
    ])

    # Define the state vector |psi> as given in the question
    # Note: 1/sqrt(2) is written as sympy.sqrt(2)/2 for precision
    psi_ket = sympy.Matrix([
        [-1/2],
        [sympy.sqrt(2)/2],
        [-1/2]
    ])

    # The multiple-choice options provided in the question
    options = {
        'A': sympy.sqrt(2) * hbar,
        'B': hbar / 2,
        'C': hbar,
        'D': hbar / sympy.sqrt(2)
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer_key = 'D'
    llm_answer_value = options[llm_answer_key]

    # --- Calculation as per Quantum Mechanics principles ---
    # The uncertainty Delta P_z is given by sqrt(<P_z^2> - <P_z>^2)

    # 1. Get the bra vector <psi|, which is the conjugate transpose of the ket vector |psi>
    psi_bra = psi_ket.H

    # 2. Verify the normalization constraint: <psi|psi> must equal 1
    norm_squared = (psi_bra * psi_ket)[0]
    if not sympy.simplify(norm_squared) == 1:
        return f"Constraint not satisfied: The state vector |psi> is not normalized. <psi|psi> = {sympy.simplify(norm_squared)}."

    # 3. Calculate the expectation value <P_z> = <psi|P_z|psi>
    exp_Pz = (psi_bra * Pz * psi_ket)[0]
    
    # 4. Calculate the operator P_z^2 by matrix multiplication
    Pz_squared = Pz * Pz

    # 5. Calculate the expectation value <P_z^2> = <psi|P_z^2|psi>
    exp_Pz_squared = (psi_bra * Pz_squared * psi_ket)[0]

    # 6. Calculate the variance, which is (<P_z^2> - <P_z>^2)
    variance = exp_Pz_squared - exp_Pz**2

    # 7. Calculate the uncertainty (standard deviation), which is the square root of the variance
    calculated_uncertainty = sympy.sqrt(variance)

    # --- Verification ---
    # Compare the independently calculated result with the LLM's answer
    if sympy.simplify(calculated_uncertainty - llm_answer_value) == 0:
        return "Correct"
    else:
        # If there is a mismatch, provide a detailed reason for the error
        reason = (
            f"Incorrect: The calculated uncertainty does not match the provided answer.\n"
            f"Provided answer ('{llm_answer_key}'): {llm_answer_value}\n"
            f"Calculated uncertainty: {calculated_uncertainty}\n"
            f"--- Intermediate Calculations ---\n"
            f"<P_z> = {exp_Pz}\n"
            f"<P_z^2> = {exp_Pz_squared}\n"
            f"Variance (<P_z^2> - <P_z>^2) = {variance}"
        )
        return reason

# Execute the check and print the result
result = check_correctness_of_answer()
print(result)