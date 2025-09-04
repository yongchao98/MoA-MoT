import numpy as np
import sympy

def check_answer():
    """
    This function checks the correctness of the provided answer for the quantum mechanics problem.
    It calculates the uncertainty Delta P_z from scratch and compares it to the given answer.
    """
    # Use sympy for symbolic representation of hbar to maintain precision.
    hbar = sympy.Symbol('hbar')

    # 1. Define the operator P_z and the state vector |psi> from the problem description.
    # The matrix form of the operator P_z.
    P_z = np.array([
        [hbar, 0, 0],
        [0, 0, 0],
        [0, 0, -hbar]
    ], dtype=object)

    # The state of the system |psi> as a column vector.
    psi_ket = np.array([
        [-1/2],
        [1/np.sqrt(2)],
        [-1/2]
    ], dtype=object)

    # 2. Verify that the state vector is normalized, i.e., <psi|psi> = 1.
    # The bra vector <psi| is the conjugate transpose of the ket vector.
    psi_bra = psi_ket.conj().T
    norm_squared = (psi_bra @ psi_ket)[0, 0]
    if not np.isclose(float(norm_squared), 1.0):
        return f"Constraint not satisfied: The state vector |psi> is not normalized. <psi|psi> = {norm_squared}, but it should be 1."

    # 3. Calculate the expectation value of P_z, <P_z> = <psi|P_z|psi>.
    Pz_psi = P_z @ psi_ket
    exp_Pz = (psi_bra @ Pz_psi)[0, 0]

    # Let's verify the intermediate calculation from the answer.
    # The answer states <P_z> = 0.
    if sympy.simplify(exp_Pz) != 0:
        return f"Incorrect intermediate calculation: The expectation value <P_z> was calculated as {sympy.simplify(exp_Pz)}, but the correct value is 0."

    # 4. Calculate the expectation value of P_z^2, <P_z^2> = <psi|P_z^2|psi>.
    # First, find the matrix for the operator P_z^2.
    Pz_squared = P_z @ P_z
    
    # The answer states P_z^2 should be diag(hbar^2, 0, hbar^2).
    expected_Pz_squared = np.array([
        [hbar**2, 0, 0],
        [0, 0, 0],
        [0, 0, hbar**2]
    ], dtype=object)
    if not np.array_equal(Pz_squared, expected_Pz_squared):
         return f"Incorrect intermediate calculation: The matrix for P_z^2 is incorrect."

    # Now calculate the expectation value.
    Pz_squared_psi = Pz_squared @ psi_ket
    exp_Pz_squared = (psi_bra @ Pz_squared_psi)[0, 0]

    # Let's verify the intermediate calculation from the answer.
    # The answer states <P_z^2> = hbar^2 / 2.
    expected_exp_Pz_squared = hbar**2 / 2
    if sympy.simplify(exp_Pz_squared - expected_exp_Pz_squared) != 0:
        return f"Incorrect intermediate calculation: The expectation value <P_z^2> was calculated as {sympy.simplify(exp_Pz_squared)}, but the correct value is {expected_exp_Pz_squared}."

    # 5. Calculate the uncertainty Delta P_z using the formula sqrt(<P_z^2> - <P_z>^2).
    variance = exp_Pz_squared - exp_Pz**2
    uncertainty = sympy.sqrt(variance)

    # 6. Compare the calculated result with the final answer.
    # The final answer is B, which corresponds to hbar / sqrt(2).
    expected_uncertainty = hbar / sympy.sqrt(2)

    if sympy.simplify(uncertainty - expected_uncertainty) == 0:
        return "Correct"
    else:
        return f"Incorrect final answer: The calculated uncertainty is {uncertainty}, but the answer B corresponds to {expected_uncertainty}."

# Run the check
result = check_answer()
print(result)