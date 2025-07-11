import math

def solve_number_of_sets():
    """
    Calculates the number of sets T with |T|=m of non-empty subsets of {1,...,n}
    such that each element i in {1,...,n} is contained in an even number of subsets in T.

    The user can modify the values of n and m below.
    """
    # Please specify the positive integers n and m for the problem.
    n = 4
    m = 5

    # The problem is for positive integers n and m.
    if n <= 0 or m <= 0:
        print("Error: n and m must be positive integers.")
        return

    # Handle small values of m as base cases.
    if m == 1 or m == 2:
        print(f"For n={n}, m={m}:")
        print("The number of such sets is 0.")
        print("<<<0>>>")
        return

    # N is the number of non-empty subsets of S.
    N = 2**n - 1

    # f is a list to store the results f_i for i from 0 to m.
    f = [0] * (m + 1)
    
    # Base case for the recurrence, f_0 = 1 (empty set sums to zero).
    f[0] = 1
    # f_1 = 0 and f_2 = 0 are implicitly set by the list initialization.

    # Iteratively compute f_i up to m using the recurrence relation.
    # The recurrence is: i * f_i = C(N, i-1) - f_{i-1} - (N - i + 2) * f_{i-2}
    for i in range(3, m + 1):
        # Using math.comb for efficient calculation of binomial coefficients.
        val_binom = math.comb(N, i - 1)
        
        # Numerator of the recurrence relation
        numerator = val_binom - f[i-1] - (N - i + 2) * f[i-2]
        
        # The result must be an integer.
        f[i] = numerator // i

    result = f[m]
    
    print(f"For n={n}, m={m}:")
    print(f"Let N = 2^n - 1 = {N}.")
    print(f"Let f(k) be the number of sets for a given size k.")
    print("The recurrence relation is k*f(k) = C(N, k-1) - f(k-1) - (N-k+2)*f(k-2), with f(0)=1, f(1)=0, f(2)=0.")

    print(f"\nFinal calculation for m = {m}:")

    val_binom_final = math.comb(N, m - 1)
    val_f_m_1 = f[m-1]
    val_f_m_2 = f[m-2]
    coeff_f_m_2 = N - m + 2
    numerator_final = val_binom_final - val_f_m_1 - coeff_f_m_2 * val_f_m_2
    
    # Printing the equation with all the numbers
    print(f"{m} * f({m}) = C({N}, {m-1}) - f({m-1}) - ({N} - {m} + 2) * f({m-2})")
    print(f"{m} * f({m}) = {val_binom_final} - {val_f_m_1} - {coeff_f_m_2} * {val_f_m_2}")
    print(f"{m} * f({m}) = {val_binom_final} - {val_f_m_1} - {coeff_f_m_2 * val_f_m_2}")
    print(f"{m} * f({m}) = {numerator_final}")
    print(f"f({m}) = {numerator_final} / {m}")
    print(f"The number of sets is {result}.")
    
    print(f"<<<{result}>>>")

solve_number_of_sets()