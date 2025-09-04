import sympy
from sympy import sin, cos, I, Matrix, Symbol, simplify, trigsimp, sqrt

def check_correctness():
    """
    Checks the correctness of the answer to the quantum mechanics problem.
    
    The problem asks for the normalized eigenvector of the operator P_n
    for the eigenvalue +hbar/2.
    
    The proposed answer is A: (cos(theta/2), sin(theta/2))
    """
    
    # 1. Define symbols and operators from the problem statement
    hbar = Symbol('hbar', real=True, positive=True)
    theta = Symbol('theta', real=True)

    # Pauli matrices scaled by hbar/2
    Px = (hbar / 2) * Matrix([[0, 1], [1, 0]])
    Pz = (hbar / 2) * Matrix([[1, 0], [0, -1]])

    # 2. Construct the operator P_n for a direction n in the x-z plane
    # The direction vector is n = (sin(theta), 0, cos(theta))
    P_n = Px * sin(theta) + Pz * cos(theta)

    # The target eigenvalue
    eigenvalue = hbar / 2

    # 3. Define the candidate eigenvector from the proposed answer (Option A)
    # Option A: (cos(theta/2), sin(theta/2))
    candidate_vector = Matrix([cos(theta/2), sin(theta/2)])

    # 4. Verify the eigenvalue equation: P_n * |psi> = lambda * |psi>
    # Calculate the left-hand side (LHS)
    lhs = P_n * candidate_vector
    # Calculate the right-hand side (RHS)
    rhs = eigenvalue * candidate_vector

    # The check is successful if (LHS - RHS) simplifies to the zero vector.
    # trigsimp is a good choice for simplifying expressions with trig functions.
    if trigsimp(lhs - rhs) != Matrix([[0], [0]]):
        return (f"Incorrect: The vector from option A, {candidate_vector.T}, is not an eigenvector "
                f"for the eigenvalue +hbar/2. The eigenvalue equation is not satisfied.")

    # 5. Verify the normalization condition: <psi|psi> = 1
    # For a real vector, this is the sum of the squares of its components.
    norm_squared = candidate_vector[0]**2 + candidate_vector[1]**2
    
    if simplify(norm_squared) != 1:
        return (f"Incorrect: The eigenvector from option A, {candidate_vector.T}, is not normalized. "
                f"Its norm squared is {simplify(norm_squared)}, not 1.")

    # 6. Verify other options are incorrect
    # Options C and D contain hbar, making them dimensionally incorrect for a state vector.
    # A state vector must be dimensionless.
    option_C_vector = Matrix([sqrt(2/3)*hbar*cos(theta/2), sqrt(1/3)*hbar*sin(theta/2)])
    if option_C_vector.has(hbar):
        # This is a simple check for the presence of the hbar symbol.
        pass # Correctly identified as dimensionally flawed.
    else:
        return "Incorrect: Logic error in checking option C's dimensions."

    # Check Option B: (cos(theta), e^(i*phi)*sin(theta)). In the x-z plane, phi=0.
    option_B_vector = Matrix([cos(theta), sin(theta)])
    lhs_B = P_n * option_B_vector
    rhs_B = eigenvalue * option_B_vector
    if trigsimp(lhs_B - rhs_B) == Matrix([[0], [0]]):
        return "Incorrect: The final answer is A, but the code found that option B is also a valid eigenvector, which contradicts the derivation."

    # If all checks pass for option A and fail for others, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)