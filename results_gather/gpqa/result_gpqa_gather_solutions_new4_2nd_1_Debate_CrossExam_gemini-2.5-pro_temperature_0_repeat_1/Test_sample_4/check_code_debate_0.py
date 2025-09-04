import sympy
from sympy import sin, cos, Matrix, I, sqrt, simplify, trigsimp

def check_answer():
    """
    This function checks the correctness of the proposed answer to the quantum mechanics problem.

    The problem asks for the normalized eigenvector of the operator P_n corresponding to the eigenvalue +ħ/2.
    P_n = (ħ/2) * [[cos(θ), sin(θ)], [sin(θ), -cos(θ)]]
    The proposed answer is (cos(θ/2), sin(θ/2)).
    """
    # Define symbolic variables
    # hbar is used for Planck's constant, theta for the angle
    hbar = sympy.Symbol('hbar', real=True, positive=True)
    theta = sympy.Symbol('theta', real=True)

    # 1. Construct the operator P_n
    # P_n = P_x*sin(theta) + P_z*cos(theta)
    # P_x = (hbar/2) * [[0, 1], [1, 0]]
    # P_z = (hbar/2) * [[1, 0], [0, -1]]
    # This simplifies to P_n = (hbar/2) * [[cos(theta), sin(theta)], [sin(theta), -cos(theta)]]
    Pn = (hbar / 2) * Matrix([
        [cos(theta), sin(theta)],
        [sin(theta), -cos(theta)]
    ])

    # 2. Define the proposed eigenvector from the final answer (Option A)
    # The final answer states the correct eigenvector is (cos(theta/2), sin(theta/2))
    psi = Matrix([
        cos(theta / 2),
        sin(theta / 2)
    ])

    # 3. Define the target eigenvalue
    eigenvalue = hbar / 2

    # 4. Check the eigenvalue equation: P_n * psi = eigenvalue * psi
    # We check if (P_n * psi) - (eigenvalue * psi) simplifies to the zero vector.
    result_vector = Pn * psi - eigenvalue * psi
    
    # Use trigsimp to simplify trigonometric expressions
    simplified_result = trigsimp(result_vector)

    # The result should be a 2x1 zero matrix
    zero_vector = Matrix([0, 0])
    if simplified_result != zero_vector:
        return (f"Incorrect: The proposed answer is not an eigenvector for the eigenvalue +ħ/2. "
                f"The result of (P_n * psi - λ * psi) should be the zero vector, but it simplified to {simplified_result}.")

    # 5. Check the normalization condition: |a|^2 + |b|^2 = 1
    # For real components, this is a^2 + b^2 = 1
    norm_squared = psi[0]**2 + psi[1]**2
    simplified_norm_squared = simplify(norm_squared)

    if simplified_norm_squared != 1:
        return (f"Incorrect: The proposed eigenvector is not normalized. "
                f"The sum of the squares of its components should be 1, but it is {simplified_norm_squared}.")

    # 6. Check other constraints mentioned in the analysis
    # Dimensionality check: The eigenvector components should be dimensionless (not contain hbar).
    # Options C and D in the prompt are:
    # C) (sqrt(2/3)*hbar, sqrt(1/3)*hbar)
    # D) (sqrt(2/3)*hbar*cos(theta/2), sqrt(1/3)*hbar*sin(theta/2))
    # We can check if hbar is a free symbol in the answer vector's components.
    if hbar in psi.free_symbols:
         return ("Incorrect: The eigenvector components should be dimensionless and not contain ħ (hbar).")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)