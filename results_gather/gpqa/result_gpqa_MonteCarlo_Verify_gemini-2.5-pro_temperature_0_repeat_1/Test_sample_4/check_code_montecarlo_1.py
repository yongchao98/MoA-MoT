import sympy
from sympy import sin, cos, I, Matrix, Symbol, simplify, trigsimp

def check_correctness():
    """
    This function symbolically verifies the eigenvector for a quantum spin operator.

    It checks two conditions for the proposed answer:
    1. It must be an eigenvector of the operator P_n with eigenvalue +hbar/2.
    2. It must be normalized to 1.
    """
    # Define symbolic variables for hbar (a real, positive constant) and theta (a real angle)
    hbar = Symbol('hbar', real=True, positive=True)
    theta = Symbol('theta', real=True)

    # Define the Pauli matrices (the core of the P operators)
    sigma_x = Matrix([[0, 1], [1, 0]])
    sigma_z = Matrix([[1, 0], [0, -1]])

    # The operator P is given by its components P_x, P_y, P_z
    # P_x = (hbar / 2) * sigma_x
    # P_z = (hbar / 2) * sigma_z

    # The direction vector n lies in the x-z plane, so n = (sin(theta), 0, cos(theta))
    n_x = sin(theta)
    n_z = cos(theta)

    # Construct the operator P_n = P_x*n_x + P_y*n_y + P_z*n_z
    # Since n_y = 0, P_n = P_x*n_x + P_z*n_z
    P_n = (hbar / 2) * (sigma_x * n_x + sigma_z * n_z)
    
    # The problem states the eigenvalue is +hbar/2
    eigenvalue = hbar / 2

    # The proposed answer from option D is the eigenvector |ψ⟩ = (cos(θ/2), sin(θ/2))
    # We represent it as a column matrix (vector)
    proposed_eigenvector = Matrix([
        cos(theta / 2),
        sin(theta / 2)
    ])

    # --- CHECK 1: Eigenvalue Equation ---
    # We need to check if P_n * |ψ⟩ = eigenvalue * |ψ⟩.
    # This is equivalent to checking if (P_n * |ψ⟩) - (eigenvalue * |ψ⟩) is the zero vector.

    # Calculate the left side of the equation: P_n * |ψ⟩
    left_side = P_n * proposed_eigenvector
    
    # Calculate the right side of the equation: eigenvalue * |ψ⟩
    right_side = eigenvalue * proposed_eigenvector

    # The check passes if the simplified difference is a zero vector.
    # trigsimp is used to simplify expressions with trigonometric functions.
    if trigsimp(left_side - right_side) != Matrix([0, 0]):
        return (f"Incorrect: The proposed eigenvector from option D does not satisfy the eigenvalue equation.\n"
                f"Applying the operator P_n to the vector yields:\n{simplify(left_side)}\n"
                f"This is not equal to the eigenvalue (+hbar/2) times the vector, which is:\n{simplify(right_side)}")

    # --- CHECK 2: Normalization ---
    # We need to check if ⟨ψ|ψ⟩ = 1.
    # For a column vector v, this is calculated as v.H * v (conjugate transpose times v).
    norm_squared = (proposed_eigenvector.H * proposed_eigenvector)[0, 0]
    
    # The check passes if the simplified norm squared is 1.
    if simplify(norm_squared) != 1:
        return (f"Incorrect: The proposed eigenvector from option D is not normalized.\n"
                f"The sum of the squares of its components, |cos(θ/2)|^2 + |sin(θ/2)|^2, simplifies to {simplify(norm_squared)}, not 1.")

    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_correctness()
print(result)