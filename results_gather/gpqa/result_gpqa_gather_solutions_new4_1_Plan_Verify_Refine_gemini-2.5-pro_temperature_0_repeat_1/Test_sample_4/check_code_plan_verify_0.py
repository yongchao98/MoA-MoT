import sympy
from sympy import sin, cos, pi, I, Matrix, Symbol, simplify, trigsimp

def check_correctness_of_answer():
    """
    This function checks the correctness of the LLM's answer to the quantum mechanics problem.

    The problem asks for the normalized eigenvector of the operator P_n = P . n
    for the eigenvalue +hbar/2, where n is a unit vector in the x-z plane.

    The LLM's final answer is D, which corresponds to the eigenvector (cos(theta/2), sin(theta/2)).
    This function will verify if this eigenvector satisfies all conditions.
    """
    try:
        # 1. Define symbols and constants from the problem statement
        hbar = Symbol('hbar', real=True, positive=True)
        theta = Symbol('theta', real=True)

        # 2. Define the operator components Px and Pz
        # Py is not needed as the direction vector n has no y-component.
        Px = (hbar / 2) * Matrix([[0, 1], [1, 0]])
        Pz = (hbar / 2) * Matrix([[1, 0], [0, -1]])

        # 3. Define the direction vector n in the x-z plane
        # n = (sin(theta), 0, cos(theta))
        nx = sin(theta)
        nz = cos(theta)

        # 4. Construct the operator P_n = P . n
        P_n = Px * nx + Pz * nz

        # 5. Define the proposed eigenvector from the LLM's answer
        # The final answer is D, which corresponds to (cos(theta/2), sin(theta/2))
        proposed_eigenvector = Matrix([
            cos(theta / 2),
            sin(theta / 2)
        ])
        
        # 6. Define the target eigenvalue
        target_eigenvalue = hbar / 2

        # 7. Check if the eigenvector is normalized
        # The norm squared should be 1. Since components are real, it's the sum of squares.
        norm_squared = proposed_eigenvector[0]**2 + proposed_eigenvector[1]**2
        if simplify(norm_squared) != 1:
            return f"Incorrect: The proposed eigenvector is not normalized. The sum of the squares of its components simplifies to {simplify(norm_squared)}, not 1."

        # 8. Check if the eigenvalue equation (P_n * psi = lambda * psi) is satisfied
        # Calculate the left-hand side (LHS): P_n * psi
        lhs = P_n * proposed_eigenvector
        # Simplify the LHS using trigonometric identities
        lhs_simplified = trigsimp(lhs)

        # Calculate the right-hand side (RHS): lambda * psi
        rhs = target_eigenvalue * proposed_eigenvector
        
        # The difference between LHS and RHS should be a zero vector
        difference = lhs_simplified - rhs
        if simplify(difference) != Matrix([[0], [0]]):
            return (f"Incorrect: The eigenvalue equation is not satisfied.\n"
                    f"Applying the operator P_n to the proposed eigenvector results in:\n{lhs_simplified}\n"
                    f"But the target (lambda * eigenvector) is:\n{rhs}")

        # 9. Check for other constraints
        # The components of an eigenvector (state vector) must be dimensionless.
        # The proposed answer (cos(theta/2), sin(theta/2)) is dimensionless.
        # Some options in the prompt contained hbar, which is dimensionally incorrect.
        # The final answer D is consistent with this principle.

        return "Correct"
    except Exception as e:
        return f"An error occurred during verification: {e}"

# Run the check
result = check_correctness_of_answer()
print(result)