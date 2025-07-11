import math

def construct_rigid_matrix(N):
    """
    This function calculates the rank parameter 'r' for a rigid matrix
    that can be constructed for a given size N using an FNP algorithm.

    The construction is theoretical and relies on a polynomial-time algorithm
    with access to an NP oracle (like a SAT solver). The core idea is:
    1. Use the NP oracle to find the truth table of a boolean function that
       is provably "hard" for small circuits. This is possible under the
       standard complexity assumption that NEXP is not in P/poly.
    2. Reshape this truth table (a long string of bits) into an N x N matrix.
    3. Theoretical results show this "hard" matrix is also rigid.

    The largest rank 'r' for which this is known to work is of the form
    r = N / polylog(N). We use r = N / (log2(N))^2 as a concrete example.
    A constant delta, e.g., delta = 0.01, is assumed.

    Args:
        N (int): The size of the N x N matrix.

    Returns:
        None. Prints the result.
    """
    if N <= 1:
        print(f"For N = {N}, the concept of rigidity is not well-defined in this context.")
        return

    # For the formula to be meaningful, N should be large enough so log(N) > 1.
    if math.log2(N) <= 1:
        print(f"N={N} is too small for this asymptotic formula.")
        return

    log_n = math.log2(N)
    # The rank 'r' is N divided by a polynomial of log(N).
    # We use (logN)^2 as a representative choice from the literature.
    r_val = N / (log_n ** 2)

    # The problem states delta is a small constant.
    delta = 0.01

    print(f"For an N x N matrix with N = {N}:")
    print(f"Given a small constant delta = {delta}, the largest rank 'r' for which we can construct a")
    print(f"(delta, r)-rigid matrix with an FNP algorithm is on the order of N/polylog(N).")
    print(f"Using the formula r = N / (log2(N))^2, we get:")
    # Outputting the numbers in the final equation as requested
    print(f"r = {N} / ({log_n:.2f})^2 = {r_val:.2f}")

# Example execution for a user-provided N.
# You can change this value to see the result for a different N.
N_input = 1024
construct_rigid_matrix(N_input)
