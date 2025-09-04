import sympy as sp

def check_eigenvector_answer():
    """
    Symbolically verifies the eigenvector for the given quantum mechanics problem.
    """
    # 1. Define symbolic variables
    # h represents h-bar (hbar), t represents theta
    h = sp.Symbol('hbar', real=True, positive=True)
    t = sp.Symbol('theta', real=True)

    # 2. Define the operator P_n in the x-z plane
    # P_n = P_x*sin(theta) + P_z*cos(theta)
    # P_x = (h/2) * [[0, 1], [1, 0]]
    # P_z = (h/2) * [[1, 0], [0, -1]]
    P_n = (h / 2) * sp.Matrix([
        [sp.cos(t), sp.sin(t)],
        [sp.sin(t), -sp.cos(t)]
    ])

    # The eigenvalue we are looking for
    eigenvalue = h / 2

    # The proposed answer from the LLM is D
    # Option D: (cos(theta/2), sin(theta/2))
    proposed_eigenvector = sp.Matrix([sp.cos(t / 2), sp.sin(t / 2)])

    # 3. Verify the eigenvalue equation: P_n * |psi> = lambda * |psi>
    
    # Calculate the left-hand side (LHS): P_n * |psi>
    lhs = P_n * proposed_eigenvector
    # Use trigonometric simplification
    lhs_simplified = sp.trigsimp(lhs)

    # Calculate the right-hand side (RHS): lambda * |psi>
    rhs = eigenvalue * proposed_eigenvector

    # Check if the equation holds by simplifying the difference
    is_eigenvector = sp.simplify(lhs_simplified - rhs) == sp.zeros(2, 1)

    # 4. Verify the normalization condition: <psi|psi> = 1
    # For a real vector, this is the dot product with itself.
    # For a complex vector, it's the conjugate transpose times the vector.
    # Here, the vector is real.
    norm_squared = proposed_eigenvector.T * proposed_eigenvector
    # The result is a 1x1 matrix, so we extract the element
    norm_squared_simplified = sp.trigsimp(norm_squared[0, 0])
    
    is_normalized = (norm_squared_simplified == 1)

    # 5. Final conclusion
    if is_eigenvector and is_normalized:
        return "Correct"
    else:
        reasons = []
        if not is_eigenvector:
            reasons.append(f"The proposed vector {proposed_eigenvector.T} is not an eigenvector for the eigenvalue +hbar/2.")
            reasons.append(f"Applying the operator P_n yields {lhs_simplified.T}, but the expected result is {rhs.T}.")
        if not is_normalized:
            reasons.append(f"The proposed vector is not normalized. Its norm squared is {norm_squared_simplified}, not 1.")
        
        # Add analysis of other options for completeness
        reasons.append("\nAnalysis of other options:")
        reasons.append("A) (cos(theta), sin(theta)) is not an eigenvector for a general theta.")
        reasons.append("B) and C) are dimensionally incorrect. A state vector's components are dimensionless probability amplitudes, but these options include hbar, giving them units of action. Their norm-squared would be proportional to hbar^2, not 1.")
        
        return "Incorrect. " + " ".join(reasons)

# Run the check and print the result
result = check_eigenvector_answer()
print(result)