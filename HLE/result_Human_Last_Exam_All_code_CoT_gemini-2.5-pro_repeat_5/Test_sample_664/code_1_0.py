import math

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def num_partitions_into_pairs(n):
    """Calculates the number of ways to partition a set of n elements into n/2 pairs."""
    if n % 2 != 0:
        return 0
    k = n // 2
    return math.factorial(n) // (math.factorial(k) * (2**k))

def calculate_involutions(n):
    """Calculates the number of involutions on a set of n elements."""
    total_involutions = 0
    print(f"Calculating the number of involutions on {n} elements (a_{n}):")
    # An involution is a permutation consisting only of 1-cycles (fixed points)
    # and 2-cycles (transpositions).
    # We sum over the number of 2-cycles, k.
    for k in range(n // 2 + 1):
        num_2_cycles = k
        num_1_cycles = n - 2 * k
        
        # Choose 2k elements for the 2-cycles out of n
        ways_to_choose = combinations(n, 2 * k)
        
        # Partition these 2k elements into k pairs
        ways_to_partition = num_partitions_into_pairs(2 * k)
        
        term = ways_to_choose * ways_to_partition
        total_involutions += term
        
        print(f"  - With {num_2_cycles} 2-cycles and {num_1_cycles} 1-cycles: C({n}, {2*k}) * {ways_to_partition} = {term}")

    return total_involutions

def calculate_both_symmetries(n):
    """
    Calculates the number of permutations symmetric to both diagonals on an n x n board.
    This corresponds to counting involutions that commute with s=(1 n)...(n/2 n/2+1).
    n must be an even number.
    """
    if n % 2 != 0:
        return 0
    
    num_s_cycles = n // 2 # s is the anti-diagonal symmetry permutation, it has n/2 2-cycles.
    total_configs = 0
    
    print(f"\nCalculating the number of configurations symmetric to BOTH diagonals:")
    print("This is the number of involutions p that commute with s=(1 8)(2 7)(3 6)(4 5).")
    print("We count based on how p permutes the 2-cycles of s.")
    
    # Case 1: 0 pairs of s-cycles are swapped by p (all 4 are fixed)
    num_swapped_pairs = 0
    num_fixed_s_cycles = num_s_cycles - 2 * num_swapped_pairs
    # For each fixed s-cycle, p has 2 choices (identity or the 2-cycle itself)
    ways_case1 = 2**num_fixed_s_cycles
    total_configs += ways_case1
    print(f"  - Case 1: 0 pairs of s-cycles swapped, {num_fixed_s_cycles} fixed.")
    print(f"    Ways = 2^{num_fixed_s_cycles} = {ways_case1}")

    # Case 2: 1 pair of s-cycles is swapped by p
    num_swapped_pairs = 1
    num_fixed_s_cycles = num_s_cycles - 2 * num_swapped_pairs
    # Choose the pair to swap
    ways_to_choose_pair = combinations(num_s_cycles, 2)
    # For the swapped pair, there are 2 ways to define p.
    # For each fixed s-cycle, there are 2 choices.
    ways_case2 = ways_to_choose_pair * (2**num_swapped_pairs) * (2**num_fixed_s_cycles)
    total_configs += ways_case2
    print(f"  - Case 2: 1 pair of s-cycles swapped, {num_fixed_s_cycles} fixed.")
    print(f"    Ways = C({num_s_cycles}, 2) * 2^{num_swapped_pairs} * 2^{num_fixed_s_cycles} = {ways_to_choose_pair} * {2**num_swapped_pairs} * {2**num_fixed_s_cycles} = {ways_case2}")

    # Case 3: 2 pairs of s-cycles are swapped by p
    num_swapped_pairs = 2
    num_fixed_s_cycles = num_s_cycles - 2 * num_swapped_pairs
    # Choose 4 cycles and partition them into 2 pairs
    ways_to_choose_pairs = num_partitions_into_pairs(num_s_cycles)
    # For each swapped pair, there are 2 ways to define p
    ways_case3 = ways_to_choose_pairs * (2**num_swapped_pairs)
    total_configs += ways_case3
    print(f"  - Case 3: 2 pairs of s-cycles swapped, {num_fixed_s_cycles} fixed.")
    print(f"    Ways = {ways_to_choose_pairs} * 2^{num_swapped_pairs} = {ways_to_choose_pairs} * {2**num_swapped_pairs} = {ways_case3}")

    return total_configs

def main():
    n = 8
    print("Step 1: Calculate configurations symmetric along the main diagonal (N_main).")
    print("This is the number of involutions on 8 elements.")
    n_main = calculate_involutions(n)
    print(f"Total number of involutions on {n} elements, N_main = {n_main}\n")
    
    print("Step 2: Calculate configurations symmetric along the anti-diagonal (N_anti).")
    print("The number of such configurations is also the number of involutions on 8 elements.")
    n_anti = n_main
    print(f"N_anti = {n_anti}\n")
    
    print("Step 3: Calculate configurations symmetric along BOTH diagonals (N_both).")
    n_both = calculate_both_symmetries(n)
    print(f"Total number of configurations symmetric to both diagonals, N_both = {n_both}\n")

    print("Step 4: Use the Principle of Inclusion-Exclusion for the final answer.")
    total_configurations = n_main + n_anti - n_both
    print(f"Total configurations = N_main + N_anti - N_both")
    print(f"Total configurations = {n_main} + {n_anti} - {n_both} = {total_configurations}")
    
    print(f"\n<<< {total_configurations} >>>")

if __name__ == "__main__":
    main()