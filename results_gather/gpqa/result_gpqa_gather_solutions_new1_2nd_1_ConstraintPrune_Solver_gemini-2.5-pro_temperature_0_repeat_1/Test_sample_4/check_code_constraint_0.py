import sympy
from sympy import sin, cos, I, Matrix, Symbol, trigsimp, simplify

def check_correctness():
    """
    Checks the correctness of the proposed answer to the quantum mechanics problem.

    The function will:
    1. Define the quantum operators and vectors using sympy.
    2. Construct the operator P_n for a direction in the x-z plane.
    3. Define the proposed eigenvector from the given answer 'C'.
    4. Verify that the eigenvector is normalized.
    5. Verify that the eigenvector satisfies the eigenvalue equation P_n * psi = (+hbar/2) * psi.
    """
    try:
        # Define symbolic variables
        hbar = Symbol('hbar', real=True, positive=True)
        theta = Symbol('theta', real=True)

        # Define the P operators based on Pauli matrices
        P_x = (hbar / 2) * Matrix([[0, 1], [1, 0]])
        P_y = (hbar / 2) * Matrix([[0, -I], [I, 0]])
        P_z = (hbar / 2) * Matrix([[1, 0], [0, -1]])

        # Define the direction vector n in the x-z plane
        n_x = sin(theta)
        n_y = 0
        n_z = cos(theta)

        # Construct the operator P_n = P.n
        P_n = P_x * n_x + P_y * n_y + P_z * n_z

        # The eigenvalue we are looking for
        eigenvalue = hbar / 2

        # The proposed answer is C: (cos(theta/2), sin(theta/2))
        # Let's check the constraints for this answer.
        # Constraint 1: Dimensionality. The components cos and sin are dimensionless. This is correct.
        # Options B and D contain hbar and are dimensionally incorrect.

        proposed_eigenvector = Matrix([cos(theta/2), sin(theta/2)])

        # Constraint 2: Normalization. The norm squared must be 1.
        # For a vector with real components, norm^2 is the sum of squares.
        norm_sq = proposed_eigenvector[0]**2 + proposed_eigenvector[1]**2
        simplified_norm_sq = simplify(norm_sq)
        if simplified_norm_sq != 1:
            return f"Incorrect: The proposed eigenvector {proposed_eigenvector.T} is not normalized. Its norm squared simplifies to {simplified_norm_sq}, not 1."

        # Constraint 3: Eigenvalue Equation. P_n * psi must equal eigenvalue * psi.
        # Let's calculate the difference and see if it's the zero vector.
        lhs = P_n * proposed_eigenvector
        rhs = eigenvalue * proposed_eigenvector
        
        difference = lhs - rhs

        # Use trigsimp to simplify trigonometric expressions
        simplified_difference = trigsimp(difference)

        if simplified_difference != Matrix([0, 0]):
            return f"Incorrect: The proposed eigenvector {proposed_eigenvector.T} does not satisfy the eigenvalue equation P_n * psi = (hbar/2) * psi. The expression (P_n*psi - (hbar/2)*psi) simplifies to {simplified_difference.T}, not the zero vector."

        # If all checks pass, the answer is correct.
        return "Correct"

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_correctness()
print(result)