import sympy

def check_correctness():
    """
    This function checks the correctness of the proposed answer to the quantum mechanics problem.
    It verifies if the given vector is indeed the normalized eigenvector of the operator P_n
    with the eigenvalue +hbar/2.
    """
    # 1. Define symbolic variables
    h_bar = sympy.Symbol('hbar', real=True, positive=True)
    theta = sympy.Symbol('theta', real=True)

    # 2. Define the operator components Px and Pz from the problem description
    # Px = (hbar/2) * Pauli_X
    Px = (h_bar / 2) * sympy.Matrix([[0, 1], [1, 0]])
    # Pz = (hbar/2) * Pauli_Z
    Pz = (h_bar / 2) * sympy.Matrix([[1, 0], [0, -1]])

    # 3. Define the unit vector n in the x-z plane
    # n = (sin(theta), 0, cos(theta)), where theta is the angle from the z-axis.
    nx = sympy.sin(theta)
    nz = sympy.cos(theta)

    # 4. Construct the operator P_n = Px*nx + Pz*nz
    P_n = Px * nx + Pz * nz

    # 5. Define the proposed eigenvector from answer B: (cos(theta/2), sin(theta/2))
    # This is represented as a column matrix (ket vector).
    proposed_eigenvector = sympy.Matrix([
        [sympy.cos(theta / 2)],
        [sympy.sin(theta / 2)]
    ])

    # 6a. Verify the Eigenvector Condition: P_n |psi> = lambda |psi>
    # The target eigenvalue is lambda = +hbar/2
    target_eigenvalue = h_bar / 2

    # Apply the operator to the proposed eigenvector
    applied_vector = P_n * proposed_eigenvector

    # The expected result from eigenvalue multiplication
    expected_vector = target_eigenvalue * proposed_eigenvector

    # Check if the two results are identical. We use simplify to handle trigonometric identities.
    is_eigenvector = sympy.simplify(applied_vector - expected_vector) == sympy.zeros(2, 1)

    if not is_eigenvector:
        return (f"Incorrect. The proposed vector is not an eigenvector for the eigenvalue +hbar/2.\n"
                f"Applying the operator P_n to the vector yields:\n{sympy.simplify(applied_vector)}\n"
                f"But the expected result (eigenvalue * vector) is:\n{sympy.simplify(expected_vector)}")

    # 6b. Verify the Normalization Condition: <psi|psi> = 1
    # The norm squared is the dot product of the conjugate transpose with the vector.
    # For a vector with real components, this is just the sum of the squares of the components.
    norm_squared = proposed_eigenvector.H * proposed_eigenvector
    
    # The result is a 1x1 matrix, so we extract the single element.
    norm_squared_val = norm_squared[0, 0]

    # Simplify the expression for the norm squared.
    is_normalized = sympy.simplify(norm_squared_val) == 1

    if not is_normalized:
        return (f"Incorrect. The proposed vector is an eigenvector with the correct eigenvalue, "
                f"but it is not normalized. Its norm squared is {sympy.simplify(norm_squared_val)}, not 1.")

    # If both conditions are met, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_correctness()
print(result)