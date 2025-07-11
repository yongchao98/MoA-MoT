import math

def solve_k_matching_complexity():
    """
    Determines the maximum integer k for which k-matchings can be counted in subcubic time
    based on the state-of-the-art algorithm and fine-grained complexity assumptions.

    The runtime complexity of counting k-matchings is given by O(n^exponent), where
    exponent = omega * floor(k/3) + (k mod 3).
    The computation is considered subcubic if exponent < 3.

    This script iterates on k to find the maximum k that satisfies this condition.
    """
    
    # The matrix multiplication exponent. The current best bound is < 2.372859 (as of 2022).
    # We use a value slightly higher but still below 2.373 for our calculations.
    omega = 2.372859
    
    print("Finding the maximum integer k for which counting k-matchings is subcubic.")
    print(f"The condition for subcubic complexity is: omega * floor(k/3) + (k mod 3) < 3")
    print(f"Using omega = {omega}\n")

    max_k_subcubic = 0
    k = 1
    
    while True:
        floor_k_div_3 = math.floor(k / 3)
        k_mod_3 = k % 3
        
        exponent = omega * floor_k_div_3 + k_mod_3
        
        is_subcubic = exponent < 3
        
        print(f"--- Checking k = {k} ---")
        print(f"Equation: {omega:.6f} * floor({k}/3) + ({k} mod 3) < 3")
        print(f"Calculation: {omega:.6f} * {floor_k_div_3} + {k_mod_3} = {exponent:.6f}")
        print(f"Result: {exponent:.6f} < 3 is {is_subcubic}")
        
        if is_subcubic:
            max_k_subcubic = k
            k += 1
            print(f"k={max_k_subcubic} is subcubic. Continuing...\n")
        else:
            print(f"\nCondition fails for k = {k}. The complexity is not subcubic.")
            print("This algorithmic bound aligns with hardness reductions, which suggest that")
            print("for k=4, the problem is not subcubic, conditioned on the hardness of counting 4-cliques.")
            break
            
    print(f"\nThe maximum integer k for which counting k-matchings is subcubic is: {max_k_subcubic}")

if __name__ == '__main__':
    solve_k_matching_complexity()
