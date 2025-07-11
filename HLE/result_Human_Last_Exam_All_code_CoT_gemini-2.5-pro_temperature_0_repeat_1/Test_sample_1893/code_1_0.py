import sympy

def calculate_neutralino_eigenvalue():
    """
    This script calculates a specific eigenvalue of the neutralino mass matrix
    under the condition of dynamic enhancement of radiative decay.

    The plan is as follows:
    1. The physical condition of "dynamic enhancement" implies that the states
       -i*gamma-tilde and H_b-tilde are mass eigenstates. This imposes the
       mathematical constraints M_1 = M_2 and cos(2*beta) = 0 on the matrix.
    2. The question asks for the eigenvalue not proportional to the adjustable
       parameters M_1, M_2, or mu. This points to a special case where these
       parameters are zero, leaving an eigenvalue that depends only on the
       constant M_Z.
    3. We construct the matrix for this specific case (M_1 = M_2 = mu = 0)
       and find its eigenvalues.
    4. The characteristic equation for the relevant part of the matrix and its
       solution will be printed.
    """
    # Define symbolic variables
    M_Z, lmbda = sympy.symbols('M_Z lambda')

    # Under the conditions M_1 = M_2 = mu = 0, the full 4x4 mass matrix simplifies.
    # The two trivial eigenvalues are 0. The non-trivial part is the 2x2
    # submatrix mixing the Zino (-i*Z-tilde) and Higgsino (H_a-tilde).
    sub_matrix = sympy.Matrix([
        [0, M_Z],
        [M_Z, 0]
    ])

    # The characteristic equation is det(sub_matrix - lambda*I) = 0
    char_eq = sub_matrix.charpoly(lmbda).as_expr()

    # The equation is lambda**2 - M_Z**2 = 0.
    # We extract the coefficients to display the numbers in the equation.
    c2 = char_eq.coeff(lmbda, 2)  # Coefficient of lambda**2
    c1 = char_eq.coeff(lmbda, 1)  # Coefficient of lambda**1
    c0_symbolic = char_eq.coeff(lmbda, 0) # Constant term is -M_Z**2
    c0_coeff = -1

    print("The characteristic equation for the relevant eigenvalues is:")
    print(f"({c2})*lambda**2 + ({c1})*lambda + ({c0_coeff})*M_Z**2 = 0")

    # Solve for the eigenvalues
    eigenvalues = sympy.solve(char_eq, lmbda)

    # The eigenvalues are M_Z and -M_Z. The question asks for "the" eigenvalue.
    # By convention for mass, we can refer to the positive value.
    positive_eigenvalue = M_Z

    print(f"\nThe non-zero eigenvalues are {eigenvalues[0]} and {eigenvalues[1]}.")
    print(f"The eigenvalue, by convention taking the positive value, is {positive_eigenvalue}.")

    # The problem implies a numerical answer. We use the known value for the Z boson mass.
    M_Z_value = 91.1876  # in GeV/c^2
    print(f"\nUsing the experimental value for the Z boson mass, M_Z ≈ {M_Z_value} GeV/c^2.")
    print(f"The computed eigenvalue is: {M_Z_value}")

if __name__ == '__main__':
    calculate_neutralino_eigenvalue()