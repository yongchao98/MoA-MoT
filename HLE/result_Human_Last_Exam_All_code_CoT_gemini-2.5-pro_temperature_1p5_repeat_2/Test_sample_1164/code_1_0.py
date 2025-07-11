def solve():
    """
    This function finds the smallest integer n >= 2 satisfying the problem's conditions
    by analyzing the different cases for the prime factors of n and n-1.
    """
    
    candidates = []

    # Case B: n-1 is a multiple of 10^9, but not in a way that violates condition 2.
    # The conditions imply n-1 = c * 10^9, where c is not a multiple of 10.
    # The smallest positive c is 1.
    n_B = 10**9 + 1
    candidates.append(n_B)
    print(f"Case B: n = 10**9 + 1")
    print(f"Smallest n from Case B is: {n_B}")
    
    # Case C: v_2(n-1)=9 and v_5(n)>=1
    # This means:
    # 1. n = 1 (mod 2**9)  => n = 1 (mod 512)
    # 2. n != 1 (mod 2**10) => n != 1 (mod 1024)
    # 3. n = 0 (mod 5)
    # From n = 5k and 5k = 1 (mod 512), we find k.
    # k = pow(5, -1, 512) = 205.
    # So k = 205 + 512j.
    # n = 5 * (205 + 512j) = 1025 + 2560j.
    # We need n != 1 (mod 1024).
    # n-1 = 1024 + 2560j = 1024 + (2 * 1024 + 512)j
    # For j=0, n=1025. n-1=1024. This means v_2(n-1)=10, which violates n!=1 (mod 1024).
    # We need v_2(n-1)=9. n-1 = 1024 + 2560j = 2**9 * (2 + 5j).
    # For v_2(n-1) to be 9, (2 + 5j) must be odd, which means j must be odd.
    # Smallest non-negative odd j is 1.
    j = 1
    n_C = 1025 + 2560 * j
    candidates.append(n_C)
    print(f"\nCase C: n = 1025 + 2560 * j, where j is the smallest non-negative odd integer.")
    print(f"j = {j}")
    print(f"n = 1025 + 2560 * {j} = {n_C}")
    print(f"Smallest n from Case C is: {n_C}")

    # Case D: v_5(n-1)=9 and v_2(n)>=1
    # This means:
    # 1. n = 1 (mod 5**9)
    # 2. n != 1 (mod 5**10)
    # 3. n = 0 (mod 2)
    # From n = 2k and 2k = 1 (mod 5**9), we find k.
    mod_5_9 = 5**9
    # k = pow(2, -1, mod_5_9) = (mod_5_9 + 1) // 2
    k = (mod_5_9 + 1) // 2
    # So k = k_sol + mod_5_9 * j
    # n = 2 * (k + mod_5_9 * j)
    # We need n-1 to not be a multiple of 5**10.
    # n-1 = 2k - 1 + 2 * mod_5_9 * j = mod_5_9 + 2 * mod_5_9 * j = (5**9)*(1+2j)
    # For v_5(n-1) to be 9, (1+2j) must not be a multiple of 5.
    # Smallest non-negative j for which this holds is j=0.
    j = 0
    n_D = 2 * (k + mod_5_9 * j)
    candidates.append(n_D)
    print(f"\nCase D: n = 2 * (({mod_5_9} + 1) / 2 + {mod_5_9} * j), where j is the smallest non-negative integer for which 1+2j is not a multiple of 5.")
    print(f"j = {j}")
    print(f"n = {n_D}")
    print(f"Smallest n from Case D is: {n_D}")
    
    # Find the overall smallest n
    result = min(candidates)
    print(f"\nComparing the candidates: {candidates}")
    print(f"The smallest positive integer n is {result}.")
    return result

solve()