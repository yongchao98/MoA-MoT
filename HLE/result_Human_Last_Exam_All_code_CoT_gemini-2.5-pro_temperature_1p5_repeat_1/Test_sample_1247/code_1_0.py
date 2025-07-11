import math

def calculate_av_333_3_1324():
    """
    Calculates the number of 1324-avoiding permutations of length 333
    with 3 inversions.
    """
    n = 333
    k = 3

    # Calculate I_n(3), the total number of permutations of length n with 3 inversions.
    # I_n(3) = (n-3) + (n-2)^2 + C(n-1, 3)
    term1 = n - 3
    term2 = (n - 2)**2
    # Using direct calculation for the binomial coefficient to show the numbers.
    comb_val = (n - 1) * (n - 2) * (n - 3) // 6
    
    i_n_k = term1 + term2 + comb_val

    # Calculate C_n^3(1324), the number of permutations of length n with 3 inversions
    # that contain the pattern 1324. The formula is n-4 for n>=4.
    c_n_k_sigma = n - 4

    # Calculate the final result
    av_n_k_sigma = i_n_k - c_n_k_sigma
    
    print("Step 1: Calculate the total number of permutations of length 333 with 3 inversions, I_333(3).")
    print(f"I_333(3) = (333 - 3) + (333 - 2)^2 + (332 * 331 * 330) / 6")
    print(f"I_333(3) = {term1} + {term2} + {comb_val}")
    print(f"I_333(3) = {i_n_k}\n")
    
    print("Step 2: Calculate the number of permutations of length 333 with 3 inversions containing the 1324 pattern, C_333^3(1324).")
    print(f"C_333^3(1324) = 333 - 4")
    print(f"C_333^3(1324) = {c_n_k_sigma}\n")
    
    print("Step 3: Calculate the number of 1324-avoiding permutations, av_333^3(1324).")
    print(f"av_333^3(1324) = I_333(3) - C_333^3(1324)")
    print(f"av_333^3(1324) = {i_n_k} - {c_n_k_sigma}")
    print(f"av_333^3(1324) = {av_n_k_sigma}")

calculate_av_333_3_1324()