import numpy as np

def solve_neutralino_eigenvalue():
    """
    Calculates the neutralino mass matrix eigenvalue based on the problem description.

    The problem asks for the eigenvalue that is not proportional to the supersymmetry
    parameters M1, M2, or mu, under conditions that lead to a "dynamic enhancement"
    of radiative decay. This enhancement occurs when the neutralino mass matrix
    becomes block-diagonal.

    This block-diagonalization isolates a 2x2 submatrix for the Zino-Higgsino system.
    In the limit where the adjustable parameters M1 and mu are zero, this submatrix's
    eigenvalues are determined solely by the Z boson mass, M_Z.
    """

    # The mass of the Z boson in GeV. We use the 2022 PDG value.
    M_Z = 91.1876

    # In the specified limit, the relevant 2x2 submatrix is:
    # M = [[0, M_Z],
    #      [M_Z, 0]]
    #
    # The characteristic equation is det(M - lambda*I) = 0, which gives:
    # (-lambda)(-lambda) - M_Z^2 = 0
    # lambda^2 - M_Z^2 = 0
    #
    # This is a quadratic equation of the form a*x^2 + b*x + c = 0.
    a = 1
    b = 0
    c = -M_Z**2

    # The problem requires outputting each number in the final equation.
    print("The characteristic equation for the eigenvalue (lambda) is:")
    print(f"({a}) * lambda^2 + ({b}) * lambda + ({c:.4f}) = 0")

    # The eigenvalues are the roots of this polynomial.
    eigenvalues = np.roots([a, b, c])

    # The two eigenvalues are +M_Z and -M_Z. The question asks for a single eigenvalue.
    # The corresponding physical mass is unique (|lambda| = M_Z). We will provide the
    # positive eigenvalue as the answer.
    positive_eigenvalue = max(eigenvalues)

    print("\nThe two eigenvalues that are not proportional to M1 or mu are +/- M_Z.")
    print(f"The positive eigenvalue is M_Z = {positive_eigenvalue:.4f} GeV.")

    # The final answer format is a number with one decimal place.
    print(f"\nFinal Answer: {positive_eigenvalue:.1f}")


if __name__ == "__main__":
    solve_neutralino_eigenvalue()