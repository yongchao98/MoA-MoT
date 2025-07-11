import math

def solve_set_problem(n, m):
    """
    Calculates the number of sets T of size m satisfying the given conditions.

    The problem is equivalent to finding the number of ways to choose m distinct
    non-zero vectors from the vector space F_2^n whose sum is the zero vector.

    This function solves it using a recurrence relation with dynamic programming.
    Let f(k) be the number of such sets of size k. The recurrence is:
    k * f(k) = C(N, k-1) - f(k-1) - (N - k + 2) * f(k-2)
    where N = 2^n - 1.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m <= 0:
        print("Error: n and m must be positive integers.")
        return

    # N is the number of non-empty subsets of S, which is |P(S)| - 1.
    N = (1 << n) - 1

    # We cannot choose m distinct elements if m is larger than the total number of elements.
    if m > N:
        print(f"For n={n}, m={m}:")
        print(f"It's impossible to choose m={m} distinct non-empty subsets from a total of {N}.")
        print("The number of such sets T is 0.")
        return

    # Use a list for dynamic programming to store the values of f(k)
    f = [0] * (m + 1)

    # Base case: f(0) = 1 (the empty set of vectors sums to zero)
    if m >= 0:
        f[0] = 1
    # Base case: f(1) = 0 (a single non-zero vector cannot sum to zero)
    # The list is already initialized to 0, so f[1] is correct.

    # Iteratively compute f(k) up to f(m)
    for k in range(2, m + 1):
        # Using math.comb for accurate and efficient calculation of combinations
        comb_val = math.comb(N, k - 1)
        
        # Recurrence relation:
        # k * f(k) = C(N, k-1) - f(k-1) - (N - k + 2) * f(k-2)
        numerator = comb_val - f[k - 1] - (N - k + 2) * f[k - 2]
        
        # The result of the division must be an integer due to the combinatorial nature of the problem.
        f[k] = numerator // k

    # Output the explanation and the final answer
    print(f"Given n = {n} and m = {m}:")
    print(f"The total number of non-empty subsets of S is N = 2^{n} - 1 = {N}.")
    print("Let f(k) be the number of valid sets T of size k.")
    print("The recurrence relation is: k * f(k) = C(N, k-1) - f(k-1) - (N - k + 2) * f(k-2)")
    print("\nCalculating the final answer f(m):")

    if m >= 2:
      final_comb = math.comb(N, m - 1)
      final_f_m_minus_1 = f[m - 1]
      final_f_m_minus_2 = f[m - 2]
      final_coeff = N - m + 2
      final_numerator = final_comb - final_f_m_minus_1 - final_coeff * final_f_m_minus_2

      print(f"f({m}) = (C({N}, {m-1}) - f({m-1}) - ({N} - {m} + 2) * f({m-2})) / {m}")
      print(f"f({m}) = ({final_comb} - {final_f_m_minus_1} - {final_coeff} * {final_f_m_minus_2}) / {m}")
      print(f"f({m}) = ({final_numerator}) / {m}")

    result = f[m]
    print(f"\nThe final number of sets T is: {result}")


# --- Example Usage ---
# You can change these values to solve for different n and m.
n_val = 4
m_val = 4

solve_set_problem(n_val, m_val)