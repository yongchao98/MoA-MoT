import sympy

def check_correctness():
    """
    This function checks the correctness of the proposed answer to the quantum mechanics problem.

    The problem asks for the normalized eigenvector of the operator P_n = P . n
    for the eigenvalue +hbar/2, where n is a unit vector in the x-z plane.

    The proposed answer is B: (cos(theta/2), sin(theta/2)).
    """
    try:
        # 1. Define the symbols and constants
        h_bar = sympy.Symbol('hbar', real=True, positive=True)
        theta = sympy.Symbol('theta', real=True)

        # 2. Define the operator components Px and Pz as given in the question
        # P_y is not needed since the direction vector n has n_y = 0
        Px = (h_bar / 2) * sympy.Matrix([[0, 1], [1, 0]])
        Pz = (h_bar / 2) * sympy.Matrix([[1, 0], [0, -1]])

        # 3. Define the direction vector n in the x-z plane
        # n = (sin(theta), 0, cos(theta))
        nx = sympy.sin(theta)
        nz = sympy.cos(theta)

        # 4. Construct the operator P_n = Px*nx + Pz*nz
        P_n = Px * nx + Pz * nz

        # 5. Define the proposed eigenvector from answer B and the eigenvalue
        # Answer B is (cos(theta/2), sin(theta/2))
        proposed_eigenvector = sympy.Matrix([
            sympy.cos(theta / 2),
            sympy.sin(theta / 2)
        ])
        
        # The eigenvalue given in the question
        proposed_eigenvalue = h_bar / 2

        # --- Verification Steps ---

        # Check 1: The eigenvector must be normalized.
        # The squared norm should be 1. The norm is the sqrt of the inner product <psi|psi>.
        # For a real vector, this is just the sum of the squares of its components.
        norm_squared = proposed_eigenvector.T * proposed_eigenvector
        # The result is a 1x1 matrix, so we extract the element.
        norm_squared_val = norm_squared[0]
        
        if sympy.simplify(norm_squared_val) != 1:
            return (f"Incorrect: The proposed eigenvector is not normalized. "
                    f"Its norm squared is {sympy.simplify(norm_squared_val)}, which must be 1.")

        # Check 2: The eigenvector must satisfy the eigenvalue equation: P_n * |psi> = lambda * |psi>
        # Calculate the left-hand side (LHS)
        lhs = P_n * proposed_eigenvector
        
        # Calculate the right-hand side (RHS)
        rhs = proposed_eigenvalue * proposed_eigenvector

        # The difference between LHS and RHS should be a zero vector.
        # We use sympy.simplify to resolve trigonometric identities.
        difference = sympy.simplify(lhs - rhs)
        
        if difference != sympy.Matrix([[0], [0]]):
            # Let's try a more powerful simplification just in case
            difference_trig = sympy.trigsimp(lhs - rhs)
            if difference_trig != sympy.Matrix([[0], [0]]):
                return (f"Incorrect: The eigenvalue equation is not satisfied.\n"
                        f"P_n * |psi> results in {sympy.trigsimp(lhs).T}\n"
                        f"lambda * |psi> results in {sympy.trigsimp(rhs).T}\n"
                        f"The difference is not the zero vector.")

        # Check 3: The eigenvector components must be dimensionless.
        # The expression (cos(theta/2), sin(theta/2)) is inherently dimensionless.
        # We can programmatically check that it does not contain hbar.
        if h_bar in proposed_eigenvector.free_symbols:
            return "Incorrect: The eigenvector components must be dimensionless, but they contain hbar."

        # If all checks pass, the answer is correct.
        return "Correct"

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check and print the result
result = check_correctness()
print(result)