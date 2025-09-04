import sympy

def check_answer():
    """
    Checks the correctness of the answer to the quantum mechanics problem.
    """
    # Define hbar as a symbolic variable
    hbar = sympy.Symbol('hbar', real=True, positive=True)
    
    # Define the operators as symbolic matrices
    P_x = (hbar / sympy.sqrt(2)) * sympy.Matrix([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])
    
    P_z = hbar * sympy.Matrix([
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ])
    
    # Define the state vector |psi>
    psi = sympy.Matrix([
        [-sympy.Rational(1, 2)],
        [1 / sympy.sqrt(2)],
        [-sympy.Rational(1, 2)]
    ])
    
    # --- Constraint 1: Check if the state vector is normalized ---
    # The bra vector <psi| is the conjugate transpose of the ket |psi>
    psi_dagger = psi.H
    norm_squared = (psi_dagger * psi)[0]
    
    if not sympy.simplify(norm_squared) == 1:
        return f"The answer is incorrect because the provided state vector is not normalized. <psi|psi> = {sympy.simplify(norm_squared)}, but it should be 1."

    # --- Constraint 2: Check if the state is an eigenstate of P_x with eigenvalue -hbar ---
    eigenvalue = -hbar
    Px_psi = P_x * psi
    expected_Px_psi = eigenvalue * psi
    
    if not sympy.simplify(Px_psi - expected_Px_psi) == sympy.zeros(3, 1):
        return f"The answer is incorrect because the state is not an eigenstate of P_x with eigenvalue -hbar as claimed. P_x|psi> is not equal to -hbar|psi>."

    # --- Constraint 3: Calculate the expectation values and uncertainty ---
    # Expectation value of P_z: <P_z> = <psi|P_z|psi>
    exp_Pz = (psi_dagger * P_z * psi)[0]
    
    # Expectation value of P_z^2: <P_z^2> = <psi|P_z^2|psi>
    Pz_squared = P_z * P_z
    exp_Pz_squared = (psi_dagger * Pz_squared * psi)[0]
    
    # Calculate the variance and uncertainty (Delta P_z)
    variance_Pz = exp_Pz_squared - exp_Pz**2
    uncertainty_Pz = sympy.sqrt(variance_Pz)
    
    # --- Constraint 4: Compare the result with the given answer (Option B) ---
    # Option B is hbar / sqrt(2)
    expected_answer = hbar / sympy.sqrt(2)
    
    if sympy.simplify(uncertainty_Pz - expected_answer) == 0:
        return "Correct"
    else:
        return f"The answer is incorrect. The calculated uncertainty is {uncertainty_Pz}, but the answer B corresponds to {expected_answer}."

# Run the check
result = check_answer()
print(result)