import sympy as sp

def check_quantum_eigenvector():
    """
    Checks if the proposed eigenvector is the correct normalized eigenvector
    for the given operator and eigenvalue.
    """
    # Define symbolic variables for ħ (hbar) and the angle θ (theta)
    hbar = sp.Symbol('hbar', real=True, positive=True)
    theta = sp.Symbol('theta', real=True)

    # Define the Pauli matrices (without the ħ/2 factor)
    sigma_x = sp.Matrix([[0, 1], [1, 0]])
    sigma_z = sp.Matrix([[1, 0], [0, -1]])

    # The direction vector n lies in the x-z plane, parameterized by θ
    # n = (sin(θ), 0, cos(θ))
    n_x = sp.sin(theta)
    n_z = sp.cos(theta)

    # Construct the operator P_n = (ħ/2) * (n_x * σ_x + n_z * σ_z)
    # P_n = (ħ/2) * (sin(θ)*σ_x + cos(θ)*σ_z)
    P_n = (hbar / 2) * (n_x * sigma_x + n_z * sigma_z)

    # The proposed answer is C: the eigenvector |ψ> = (cos(θ/2), sin(θ/2))
    # We represent it as a column matrix (ket vector)
    psi = sp.Matrix([
        sp.cos(theta / 2),
        sp.sin(theta / 2)
    ])

    # The given eigenvalue is λ = +ħ/2
    eigenvalue = hbar / 2

    # --- Constraint 1: Check if the eigenvector is normalized ---
    # The norm-squared is <ψ|ψ> = ψ† * ψ, where ψ† is the conjugate transpose.
    # Since the components are real, this is just ψ.T * ψ.
    norm_squared = psi.T * psi
    # The result is a 1x1 matrix, so we extract the element.
    # We use simplify to resolve the trigonometric identity cos²(x) + sin²(x) = 1.
    simplified_norm = sp.simplify(norm_squared[0, 0])

    if simplified_norm != 1:
        return (f"Incorrect: The proposed eigenvector is not normalized. "
                f"Its norm-squared is {simplified_norm}, but it should be 1.")

    # --- Constraint 2: Check if it satisfies the eigenvalue equation P_n|ψ> = λ|ψ> ---
    # Calculate the left-hand side (LHS): P_n * |ψ>
    lhs = P_n * psi
    # Calculate the right-hand side (RHS): λ * |ψ>
    rhs = eigenvalue * psi

    # Check if LHS equals RHS by simplifying their difference.
    # If the result is a zero vector, the equation holds true.
    difference = lhs - rhs
    simplified_difference = sp.simplify(difference)

    if simplified_difference != sp.Matrix([0, 0]):
        # The simplification uses trigonometric identities like:
        # sin(θ) = 2*sin(θ/2)*cos(θ/2)
        # cos(θ) = cos²(θ/2) - sin²(θ/2)
        return (f"Incorrect: The proposed eigenvector does not satisfy the eigenvalue equation. "
                f"P_n*|ψ> - λ*|ψ> should be a zero vector, but it simplifies to:\n{simplified_difference}")

    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_quantum_eigenvector()
print(result)