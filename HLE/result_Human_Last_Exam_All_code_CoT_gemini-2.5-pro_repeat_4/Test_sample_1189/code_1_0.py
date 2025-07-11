import math

def solve_sets_problem(n, m):
    """
    Calculates the number of sets T satisfying the given conditions.
    
    Args:
        n (int): The size of the base set S, a positive integer.
        m (int): The size of the set T, a positive integer.
    """
    
    # Step 1: Initialize constants and base cases for the recurrence.
    N = 2**n - 1
    
    if m < 0:
        print("Error: m must be a non-negative integer.")
        return
    
    if n <= 0:
        print("Error: n must be a positive integer.")
        return

    # Handle the case where m is larger than the number of available subsets.
    if m > N:
        print(f"For n={n}, m={m}:")
        print(f"The total number of non-empty subsets is N = 2^{n}-1 = {N}.")
        print(f"Since m > N, it is impossible to choose {m} distinct subsets. The answer is 0.")
        print("<<<0>>>")
        return

    # memo_f stores the values of f(k), our target function.
    memo_f = {0: 1, 1: 0}
    
    # Step 2: Iteratively compute f(k) up to m using the recurrence relation.
    for k in range(2, m + 1):
        # f(k) = (C(N, k-1) - f(k-1) - (N - k + 2) * f(k-2)) / k
        comb_val = math.comb(N, k - 1)
        f_k_minus_1 = memo_f[k - 1]
        f_k_minus_2 = memo_f[k - 2]
        
        numerator = comb_val - f_k_minus_1 - (N - k + 2) * f_k_minus_2
        
        # The result of the division is guaranteed to be an integer.
        memo_f[k] = numerator // k

    result = memo_f[m]

    # Step 3: Print the final calculation and result as requested.
    print(f"For n={n}, m={m}:")
    print(f"Let f(k) be the number of valid sets of size k. The total number of available non-empty subsets of S is N = 2^n - 1 = {N}.")
    
    if m >= 2:
        comb_val_final = math.comb(N, m - 1)
        f_m_minus_1 = memo_f[m - 1]
        f_m_minus_2 = memo_f[m - 2]
        N_m_2_factor = N - m + 2
        
        print("\nThe recurrence relation used is:")
        print(f"k * f(k) = C(N, k-1) - f(k-1) - (N - k + 2) * f(k-2)")
        
        print("\nThe final calculation step for f(m) is:")
        print(f"{m} * f({m}) = C({N}, {m-1}) - f({m-1}) - ({N} - {m} + 2) * f({m-2})")
        print(f"{m} * f({m}) = {comb_val_final} - {f_m_minus_1} - ({N_m_2_factor}) * {f_m_minus_2}")
        term3 = N_m_2_factor * f_m_minus_2
        print(f"{m} * f({m}) = {comb_val_final} - {f_m_minus_1} - {term3}")
        numerator_final = comb_val_final - f_m_minus_1 - term3
        print(f"{m} * f({m}) = {numerator_final}")
        print(f"f({m}) = {numerator_final} // {m} = {result}")

    elif m == 1:
        print("\nBy definition, f(1) = 0, as a single non-empty subset cannot satisfy the condition for its own elements.")
    elif m == 0:
        print("\nBy definition, f(0) = 1, representing the empty collection of subsets, which trivially satisfies the conditions.")

    print(f"<<<{result}>>>")


if __name__ == '__main__':
    # --- Set the values for n and m here ---
    n_val = 4
    m_val = 4
    # -----------------------------------------
    
    solve_sets_problem(n_val, m_val)