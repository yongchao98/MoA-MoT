import sympy

def check_correctness():
    """
    This function checks the correctness of the provided answer to the quantum mechanics problem.

    The problem asks for the normalized eigenvector of the operator P_n for the eigenvalue +hbar/2.
    The final answer provided by the LLM is C, which corresponds to the vector (cos(theta/2), sin(theta/2)).

    The checker will perform two main verifications:
    1.  Check if the proposed eigenvector is normalized (i.e., its norm is 1).
    2.  Check if the proposed eigenvector satisfies the eigenvalue equation: P_n |psi> = lambda |psi>.
    """
    try:
        # Define symbolic variables for Planck's constant and the angle
        h_bar = sympy.Symbol('hbar', real=True, positive=True)
        theta = sympy.Symbol('theta', real=True)

        # --- Step 1: Construct the operator P_n ---
        # Define the operator components Px and Pz from the problem description
        Px = (h_bar / 2) * sympy.Matrix([[0, 1], [1, 0]])
        Pz = (h_bar / 2) * sympy.Matrix([[1, 0], [0, -1]])

        # Define the components of the unit vector n in the x-z plane
        nx = sympy.sin(theta)
        nz = sympy.cos(theta)

        # Construct the operator P_n = Px*nx + Pz*nz
        P_n = Px * nx + Pz * nz

        # --- Step 2: Define the candidate eigenvector and eigenvalue ---
        # The candidate answer is C: (cos(theta/2), sin(theta/2))
        psi_C = sympy.Matrix([sympy.cos(theta/2), sympy.sin(theta/2)])
        
        # The given eigenvalue
        eigenvalue = h_bar / 2

        # --- Step 3: Verify the constraints ---

        # Constraint 1: The eigenvector must be normalized.
        # The squared norm of the vector should be 1.
        # For a vector with real components, the norm squared is the dot product with itself.
        norm_squared = psi_C.dot(psi_C)
        if sympy.simplify(norm_squared) != 1:
            return f"Incorrect. The candidate eigenvector {psi_C.T} from option C is not normalized. Its norm squared simplifies to {sympy.simplify(norm_squared)}, not 1."

        # Constraint 2: The eigenvector must satisfy the eigenvalue equation P_n |psi> = lambda |psi>.
        # This is equivalent to checking if (P_n - lambda*I) |psi> = 0.
        
        # Calculate LHS = P_n * |psi>
        lhs = P_n * psi_C
        
        # Calculate RHS = lambda * |psi>
        rhs = eigenvalue * psi_C
        
        # The difference should be a zero vector. We use trigsimp to handle the identities.
        difference = sympy.trigsimp(lhs - rhs)
        
        if difference != sympy.zeros(2, 1):
            return f"Incorrect. The candidate eigenvector {psi_C.T} from option C does not satisfy the eigenvalue equation. (P_n * psi) simplifies to {sympy.trigsimp(lhs).T}, while (lambda * psi) is {rhs.T}. The difference is not the zero vector."

        # Constraint 3: Dimensionality check of other options.
        # Options B and D contain h_bar, which has units of action.
        # State vectors in quantum mechanics are dimensionless, as their squared components represent probabilities.
        # Therefore, options B and D are incorrect on dimensional grounds.
        
        # If all checks pass, the answer is correct.
        return "Correct"

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result
result = check_correctness()
print(result)