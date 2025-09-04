import sympy

def check_correctness():
    """
    This function checks if the proposed answer B is the correct normalized eigenvector.
    """
    # Define the symbolic constants from the problem
    h_bar = sympy.Symbol('hbar', real=True, positive=True)
    theta = sympy.Symbol('theta', real=True)

    # Define the operator components Px and Pz as given in the problem
    # Px = (hbar/2) * [[0, 1], [1, 0]]
    Px = (h_bar / 2) * sympy.Matrix([[0, 1], [1, 0]])
    # Pz = (hbar/2) * [[1, 0], [0, -1]]
    Pz = (h_bar / 2) * sympy.Matrix([[1, 0], [0, -1]])

    # The direction vector n lies in the x-z plane.
    # A unit vector in this plane can be parameterized by the polar angle theta.
    # n = (sin(theta), 0, cos(theta))
    nx = sympy.sin(theta)
    nz = sympy.cos(theta)

    # Construct the operator P_n = Px*nx + Pz*nz
    P_n = Px * nx + Pz * nz

    # The proposed answer is B) (cos(theta/2), sin(theta/2))
    # We represent this as a column vector (a "ket" in quantum mechanics)
    proposed_eigenvector = sympy.Matrix([
        [sympy.cos(theta / 2)],
        [sympy.sin(theta / 2)]
    ])

    # --- CHECK 1: Eigenvector condition ---
    # We must verify that P_n * |psi> = (+hbar/2) * |psi>
    
    # Apply the operator P_n to the proposed eigenvector
    action_on_vector = P_n * proposed_eigenvector
    
    # This is what the result should be if the eigenvalue is +hbar/2
    expected_result = (h_bar / 2) * proposed_eigenvector

    # The two expressions must be identical. We can check this by simplifying their difference.
    # If the difference simplifies to a zero vector, the condition is met.
    # We use trigsimp to handle the trigonometric identities involved.
    difference = sympy.trigsimp(action_on_vector - expected_result)
    
    is_eigenvector = (difference == sympy.zeros(2, 1))

    if not is_eigenvector:
        return f"Incorrect. The proposed vector {proposed_eigenvector.T} is not an eigenvector of P_n with eigenvalue +hbar/2. The eigenvalue equation is not satisfied."

    # --- CHECK 2: Normalization condition ---
    # We must verify that the vector is normalized, i.e., its norm is 1.
    # The squared norm is the inner product of the vector with itself, which is
    # the conjugate transpose of the vector multiplied by the vector.
    # For a real vector, this is just the transpose multiplied by the vector.
    norm_squared = proposed_eigenvector.T * proposed_eigenvector
    
    # The result is a 1x1 matrix containing the value of the norm squared.
    # We simplify it to check if it equals 1.
    simplified_norm_squared = sympy.simplify(norm_squared[0, 0])

    is_normalized = (simplified_norm_squared == 1)

    if not is_normalized:
        return f"Incorrect. The proposed eigenvector {proposed_eigenvector.T} is not normalized. Its norm squared is {simplified_norm_squared}, but it should be 1."

    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_correctness()
print(result)