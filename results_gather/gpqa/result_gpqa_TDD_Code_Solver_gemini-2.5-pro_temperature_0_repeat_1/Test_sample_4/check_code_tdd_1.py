import sympy

def check_spin_eigenvector():
    """
    This function checks the correctness of the proposed eigenvector for a spin-1/2 particle.
    
    Problem details:
    - Operator: P_n = P . n, where P is the spin operator (hbar/2 * Pauli matrices).
    - Direction n: An arbitrary direction in the x-z plane, n = (sin(theta), 0, cos(theta)).
    - Eigenvalue: +hbar/2.
    - Proposed eigenvector (Answer B): |psi> = [cos(theta/2), sin(theta/2)]^T.
    
    The function performs two checks:
    1. Eigenvalue check: Verifies if P_n |psi> = (+hbar/2) |psi>.
    2. Normalization check: Verifies if <psi|psi> = 1.
    
    Returns:
        str: "Correct" if both checks pass, otherwise a string explaining the failure.
    """
    
    # Define symbolic variables and constants
    hbar = sympy.Symbol('hbar', real=True, positive=True)
    theta = sympy.Symbol('theta', real=True)
    i = sympy.I

    # Define the spin operator components (P_x, P_y, P_z)
    # These are hbar/2 times the Pauli matrices.
    P_x = (hbar / 2) * sympy.Matrix([[0, 1], [1, 0]])
    P_y = (hbar / 2) * sympy.Matrix([[0, -i], [i, 0]])
    P_z = (hbar / 2) * sympy.Matrix([[1, 0], [0, -1]])

    # Define the direction vector n in the x-z plane
    # n = (sin(theta), 0, cos(theta))
    n_x = sympy.sin(theta)
    n_y = 0
    n_z = sympy.cos(theta)

    # Construct the operator P_n = P_x*n_x + P_y*n_y + P_z*n_z
    P_n = P_x * n_x + P_y * n_y + P_z * n_z
    
    # The eigenvalue to check against is +hbar/2
    eigenvalue = hbar / 2

    # The proposed eigenvector from option B is [cos(theta/2), sin(theta/2)]
    proposed_eigenvector = sympy.Matrix([
        sympy.cos(theta / 2),
        sympy.sin(theta / 2)
    ])

    # --- CHECK 1: Eigenvalue Equation Constraint ---
    # We need to verify that P_n * |psi> = eigenvalue * |psi>
    
    # Calculate the left-hand side (LHS) of the equation
    lhs = P_n * proposed_eigenvector
    
    # Calculate the right-hand side (RHS) of the equation
    rhs = eigenvalue * proposed_eigenvector

    # Check if LHS and RHS are equal by simplifying their difference.
    # If the difference is a zero vector, the equation holds.
    # sympy.simplify() is crucial here to resolve the trigonometric identities.
    if sympy.simplify(lhs - rhs) != sympy.zeros(2, 1):
        return (f"Incorrect: The proposed vector is not an eigenvector for the eigenvalue +hbar/2. "
                f"Applying the operator P_n to the vector yields {sympy.simplify(lhs).T}, "
                f"but it should be {sympy.simplify(rhs).T}.")

    # --- CHECK 2: Normalization Constraint ---
    # The eigenvector must be normalized, i.e., <psi|psi> = 1.
    # <psi|psi> is calculated as the conjugate transpose of the vector times the vector.
    # Since the vector is real, the conjugate transpose (H) is the same as the transpose.
    norm_squared = (proposed_eigenvector.H * proposed_eigenvector)[0, 0]
    
    # Simplify the result. We expect cos(x)^2 + sin(x)^2 = 1.
    if sympy.simplify(norm_squared) != 1:
        return (f"Incorrect: The proposed eigenvector is not normalized. "
                f"Its norm squared <psi|psi> is {sympy.simplify(norm_squared)}, but it should be 1.")

    # If both checks pass, the answer is correct.
    return "Correct"

# You can run the check by calling the function and printing its output.
# For example:
# result = check_spin_eigenvector()
# print(result)
# This will print "Correct" if the provided answer B is valid.