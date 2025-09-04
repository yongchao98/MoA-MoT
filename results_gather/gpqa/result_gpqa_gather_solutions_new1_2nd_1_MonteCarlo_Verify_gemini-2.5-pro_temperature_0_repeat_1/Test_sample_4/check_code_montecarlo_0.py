import sympy

def check_answer():
    """
    Checks the correctness of the answer to the quantum mechanics problem.
    
    The problem asks for the normalized eigenvector of the operator P_n 
    corresponding to the eigenvalue +hbar/2.
    """
    # Define symbolic variables
    # hbar is a physical constant, theta is the angle
    hbar, theta = sympy.symbols('hbar theta', real=True, positive=True)
    # i is the imaginary unit
    i = sympy.I

    # Define the Pauli matrices scaled by hbar/2 as given in the question
    P_x = (hbar / 2) * sympy.Matrix([[0, 1], [1, 0]])
    P_y = (hbar / 2) * sympy.Matrix([[0, -i], [i, 0]])
    P_z = (hbar / 2) * sympy.Matrix([[1, 0], [0, -1]])

    # Define the direction vector n in the x-z plane
    # n = (sin(theta), 0, cos(theta))
    n_x = sympy.sin(theta)
    n_y = 0
    n_z = sympy.cos(theta)

    # Construct the operator P_n = P.n
    P_n = P_x * n_x + P_y * n_y + P_z * n_z

    # The eigenvalue we are interested in
    eigenvalue = hbar / 2

    # The options as provided in the question text
    # Note: The problem states n is in the x-z plane, so the azimuthal angle phi is 0.
    # This makes exp(i*phi) = 1.
    options = {
        'A': sympy.Matrix([sympy.cos(theta), sympy.sin(theta)]),
        'B': "Contains hbar, dimensionally incorrect",
        'C': sympy.Matrix([sympy.cos(theta/2), sympy.sin(theta/2)]),
        'D': "Contains hbar, dimensionally incorrect"
    }
    
    # The proposed answer is C
    proposed_answer_key = 'C'
    psi_proposed = options[proposed_answer_key]

    # --- Verification Step 1: Dimensionality Check ---
    # Options B and D contain hbar, which has units of angular momentum.
    # The components of a normalized state vector must be dimensionless.
    # This is a conceptual check.
    if proposed_answer_key in ['B', 'D']:
        return f"Incorrect. The answer '{proposed_answer_key}' is dimensionally incorrect because its components contain hbar, which has physical units. State vectors must be dimensionless."

    # --- Verification Step 2: Eigenvalue Equation Check ---
    # We need to check if P_n * psi = eigenvalue * psi
    lhs = P_n * psi_proposed
    rhs = eigenvalue * psi_proposed

    # Simplify the difference. If it's a zero vector, the equation holds.
    difference = sympy.simplify(lhs - rhs)
    
    if not difference.is_zero_matrix:
        # Check other options to see if any of them work
        for key, val in options.items():
            if isinstance(val, sympy.Matrix):
                diff_other = sympy.simplify(P_n * val - eigenvalue * val)
                if diff_other.is_zero_matrix:
                     return f"Incorrect. The proposed answer '{proposed_answer_key}' is not the correct eigenvector. The correct eigenvector is option '{key}': {val.T}."

        return f"Incorrect. The proposed answer '{proposed_answer_key}' does not satisfy the eigenvalue equation P_n|psi> = (hbar/2)|psi>. The result of P_n|psi> - (hbar/2)|psi> is {difference.T}, not the zero vector."

    # --- Verification Step 3: Normalization Check ---
    # The norm-squared should be 1. For a vector [a, b], this is a*conj(a) + b*conj(b).
    norm_squared = psi_proposed[0]**2 + psi_proposed[1]**2
    
    if sympy.simplify(norm_squared) != 1:
        return f"Incorrect. The proposed answer '{proposed_answer_key}' satisfies the eigenvalue equation but is not normalized. Its norm-squared is {sympy.simplify(norm_squared)}, not 1."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)