import math

def max_complex_zeros(N):
    """
    Calculates the maximal possible number of complex zeros for the given problem.

    A zero k_j is complex if Re(k_j) != 0 and Im(k_j) != 0.
    This condition on all k_j's is met if k_1 is a complex number
    not on the real or imaginary axis.

    Args:
        N: The dimension of the matrix, an integer >= 1.

    Returns:
        The maximal number of such complex zeros.
    """
    if not isinstance(N, int) or N < 1:
        raise ValueError("N must be an integer greater than or equal to 1.")

    # Case N=1:
    # The equation is A_11 + k_1 = 0, which has one real root k_1 = -A_11.
    # The degree of the polynomial is N * 2^(N-1) = 1 * 2^0 = 1.
    # A real polynomial of odd degree must have at least one real root.
    # So, for N=1, there is always one real root, and 0 complex roots.
    if N == 1:
        print(f"For N = 1:")
        print("Total zeros = 1 * 2^(1-1) = 1. This polynomial is linear and must have one real root.")
        result = 0
        print(f"Maximal number of complex zeros = {result}")
        print("-" * 20)
        return

    total_zeros = N * (2**(N - 1))
    
    # Case N is even:
    # At least one real root is guaranteed for each of the 2^(N-1) sign choices
    # by the Intermediate Value Theorem.
    # Assuming these roots are distinct, we have at least 2^(N-1) real roots.
    if N % 2 == 0:
        print(f"For N = {N} (even):")
        forced_real_roots = 2**(N - 1)
        result = total_zeros - forced_real_roots
        
        # Breakdown of the calculation
        # (N-1) * 2^(N-1)
        term1 = N - 1
        term2 = 2**(N - 1)

        print(f"Total zeros = {N} * 2^({N}-1) = {total_zeros}")
        print(f"Minimal number of real zeros = 2^({N}-1) = {forced_real_roots}")
        print(f"Maximal number of complex zeros = Total - Real = {total_zeros} - {forced_real_roots} = {result}")
        print(f"Which is also = ({N}-1) * 2^({N}-1) = {term1} * {term2} = {result}")

    # Case N is odd and > 1:
    # The degree N * 2^(N-1) is even. There is no theorem forcing real or purely
    # imaginary roots for a generic choice of A. Thus, all roots can be complex.
    else:
        print(f"For N = {N} (odd, >1):")
        result = total_zeros
        print(f"Total zeros = {N} * 2^({N}-1) = {result}")
        print("For a generic choice of A, all zeros can be complex.")
        print(f"Maximal number of complex zeros = {result}")

    print("-" * 20)

# Example calculations for N from 1 to 5
max_complex_zeros(1)
max_complex_zeros(2)
max_complex_zeros(3)
max_complex_zeros(4)
max_complex_zeros(5)

# Example for a larger N
max_complex_zeros(10)
