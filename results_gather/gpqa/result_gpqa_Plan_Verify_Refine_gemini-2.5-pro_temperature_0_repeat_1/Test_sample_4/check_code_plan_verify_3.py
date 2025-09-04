import sympy

def check_correctness():
    """
    This function checks if the proposed answer is the correct normalized eigenvector.
    It verifies two conditions:
    1. The vector is normalized.
    2. The vector satisfies the eigenvalue equation P_n |psi> = (+hbar/2) |psi>.
    """
    # Define the symbolic variables for Planck's constant and the angle theta.
    # hbar is real and positive. theta is a real angle.
    h_bar = sympy.Symbol('hbar', real=True, positive=True)
    theta = sympy.Symbol('theta', real=True)

    # 1. Construct the operator P_n from the problem description.
    # Define the Pauli matrices multiplied by hbar/2 as given in the question.
    Px = (h_bar / 2) * sympy.Matrix([[0, 1], [1, 0]])
    Pz = (h_bar / 2) * sympy.Matrix([[1, 0], [0, -1]])
    
    # A unit vector n in the x-z plane can be parameterized by an angle theta
    # with the z-axis: n = (sin(theta), 0, cos(theta)).
    nx = sympy.sin(theta)
    nz = sympy.cos(theta)
    
    # The operator for spin along direction n is P_n = Px*nx + Pz*nz.
    P_n = Px * nx + Pz * nz

    # 2. Define the proposed eigenvector from answer B.
    # The proposed eigenvector is |psi> = [cos(theta/2), sin(theta/2)]
    proposed_eigenvector = sympy.Matrix([
        sympy.cos(theta/2),
        sympy.sin(theta/2)
    ])

    # 3. Check if the proposed eigenvector is normalized.
    # The squared norm is the sum of the squares of the absolute values of its components.
    # For a real vector, this is v.T * v.
    norm_squared = proposed_eigenvector.T * proposed_eigenvector
    # The result is a 1x1 matrix, so we extract the element.
    simplified_norm_squared = sympy.simplify(norm_squared[0])
    
    if simplified_norm_squared != 1:
        return f"Incorrect. The proposed eigenvector is not normalized. Its norm squared is {simplified_norm_squared}, but it should be 1."

    # 4. Check if the proposed eigenvector satisfies the eigenvalue equation.
    # We need to verify if P_n * |psi> = (+hbar/2) * |psi>.
    # This is equivalent to checking if (P_n - (+hbar/2)*I) * |psi> = 0.
    
    # The eigenvalue to check against.
    eigenvalue = h_bar / 2
    
    # Calculate the action of the operator on the vector.
    action_of_P_on_psi = P_n * proposed_eigenvector
    
    # Calculate the eigenvalue times the vector.
    eigenvalue_times_psi = eigenvalue * proposed_eigenvector
    
    # The difference between the two sides of the equation should be a zero vector.
    # We use sympy.simplify to reduce the trigonometric expressions.
    difference = sympy.simplify(action_of_P_on_psi - eigenvalue_times_psi)
    
    # The zero vector in 2 dimensions.
    zero_vector = sympy.zeros(2, 1)
    
    if difference != zero_vector:
        return f"Incorrect. The proposed vector does not satisfy the eigenvalue equation. P_n*|psi> - lambda*|psi> simplifies to {difference}, not the zero vector."

    # 5. If both checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_correctness()
print(result)