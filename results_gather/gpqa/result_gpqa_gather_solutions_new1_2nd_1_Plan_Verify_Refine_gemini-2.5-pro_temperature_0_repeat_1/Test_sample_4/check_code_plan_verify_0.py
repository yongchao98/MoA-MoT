import sympy
from sympy import Matrix, symbols, cos, sin, I, sqrt, simplify, trigsimp

def check_answer():
    """
    Checks the correctness of the proposed answer for the quantum mechanics eigenvector problem.
    """
    # Define symbolic variables
    # hbar is a real, positive constant (Planck's constant over 2*pi)
    # theta is the real-valued angle of the direction vector
    hbar = symbols('hbar', real=True, positive=True)
    theta = symbols('theta', real=True)

    # Define the operator components based on the problem description
    # P_x = (hbar/2) * sigma_x
    P_x = (hbar / 2) * Matrix([[0, 1], [1, 0]])
    # P_z = (hbar/2) * sigma_z
    P_z = (hbar / 2) * Matrix([[1, 0], [0, -1]])

    # Define the direction vector n in the x-z plane
    n_x = sin(theta)
    n_z = cos(theta)

    # Construct the operator P_n = P_x*n_x + P_z*n_z
    P_n = P_x * n_x + P_z * n_z

    # The eigenvalue we are looking for
    eigenvalue = hbar / 2

    # The final answer provided is B, which corresponds to the eigenvector (cos(theta/2), sin(theta/2))
    # Let's define this as a column vector
    proposed_eigenvector = Matrix([
        cos(theta / 2),
        sin(theta / 2)
    ])

    # --- Verification Steps ---

    # 1. Check for correct dimensionality. A normalized state vector must be dimensionless.
    #    Its components should not contain hbar.
    if any(c.has(hbar) for c in proposed_eigenvector):
        return f"Incorrect: The components of a normalized eigenvector must be dimensionless, but the proposed answer contains hbar. The components are {proposed_eigenvector.T}."

    # 2. Check for normalization. The sum of the squares of the components must equal 1.
    #    For a complex vector v = [a, b], norm is a*conj(a) + b*conj(b). Here, components are real.
    norm_squared = proposed_eigenvector[0]**2 + proposed_eigenvector[1]**2
    if simplify(norm_squared) != 1:
        return f"Incorrect: The proposed eigenvector is not normalized. The sum of the squares of its components is {simplify(norm_squared)}, not 1."

    # 3. Check if it satisfies the eigenvalue equation: P_n * v = lambda * v
    #    Calculate the left-hand side (LHS)
    lhs = P_n * proposed_eigenvector
    # Calculate the right-hand side (RHS)
    rhs = eigenvalue * proposed_eigenvector

    # The difference between LHS and RHS should be a zero vector.
    # We use trigsimp to simplify trigonometric expressions.
    difference = trigsimp(lhs - rhs)

    if difference != Matrix([[0], [0]]):
        # Let's show the calculation to see where it fails
        # Simplified LHS:
        # P_n * v = hbar/2 * [[cos(t), sin(t)], [sin(t), -cos(t)]] * [cos(t/2), sin(t/2)]
        #         = hbar/2 * [cos(t)cos(t/2)+sin(t)sin(t/2), sin(t)cos(t/2)-cos(t)sin(t/2)]
        #         = hbar/2 * [cos(t-t/2), sin(t-t/2)]
        #         = hbar/2 * [cos(t/2), sin(t/2)]
        # This is exactly the RHS. The check should pass. If it fails, there's an issue.
        return f"Incorrect: The proposed eigenvector does not satisfy the eigenvalue equation P_n * v = (hbar/2) * v. The difference (LHS - RHS) is {difference.T}, not [0, 0]."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)