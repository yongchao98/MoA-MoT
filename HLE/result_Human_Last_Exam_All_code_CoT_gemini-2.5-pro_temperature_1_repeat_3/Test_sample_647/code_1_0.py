import math

def solve_committee_election():
    """
    Solves the committee election problem by calculating the satisfaction for
    core and EJR committees and then finding their ratio.
    """
    
    # 1. Define problem parameters
    num_voters_total = 10
    num_voters_N = 8
    committee_size_k = 20

    print("--- Problem Analysis ---")
    print(f"Total voters (n): {num_voters_total}, Committee size (k): {committee_size_k}, Group N size: {num_voters_N}")
    print("-" * 30)

    # 2. Calculate s(N, W_1) for the Core
    print("--- Part 1: Lowest satisfaction for a Core committee (W_1) ---")
    
    # The Core property guarantees a proportional share for any group of voters S.
    # Share for S = (|S|/n) * k
    # For group N={1..8}, share = (8/10)*20 = 16. W must have >= 16 members from {1..24}.
    # For voter {9}, share = (1/10)*20 = 2. W must have >= 2 members from {25..28}.
    # For voter {10}, share = (1/10)*20 = 2. W must have >= 2 members from {29..32}.
    
    w_N_core = 16
    w_9_core = 2
    w_10_core = 2
    
    print(f"A core committee must have at least {w_N_core} members from the set of candidates approved by group N, {{1..24}}.")
    print(f"It must also have at least {w_9_core} from voter 9's set and {w_10_core} from voter 10's set.")
    print(f"Since {w_N_core} + {w_9_core} + {w_10_core} = {committee_size_k}, this distribution is exact.")

    # Satisfaction for group N is s(N, W) = 8 * w_common + w_pairs,
    # where w_common = |W intersect {1..8}| and w_pairs = |W intersect {9..24}|.
    # We know w_common + w_pairs = w_N_core = 16.
    # So, s(N, W) = 8*w_common + (16 - w_common) = 7*w_common + 16.
    
    # To minimize satisfaction, we minimize w_common.
    # We can choose all 16 members from {9..24}, so min w_common is 0.
    w_common_min_core = 0
    s_N_W1 = 7 * w_common_min_core + 16

    print(f"To minimize satisfaction for N, we select the minimum possible number of common candidates ({w_common_min_core}).")
    print(f"s(N, W_1) = 7 * {w_common_min_core} + 16 = {s_N_W1}")
    print("-" * 30)
    
    # 3. Calculate s(N, W_2) for EJR
    print("--- Part 2: Lowest satisfaction for an EJR committee (W_2) ---")
    
    # The key EJR constraint comes from any 4-voter subgroup of N.
    # For such a group, l=4, l/n*k = 8, and they share 8 common candidates {1..8}.
    # EJR requires that for any 4 voters from N, one must have satisfaction >= 8.
    # This implies at most 3 voters in N can have satisfaction < 8.
    # So, at least 5 voters in N must have satisfaction >= 8.
    
    print("EJR requires that at least 5 voters in group N achieve a satisfaction of 8 or more.")
    
    # To minimize s(N, W) = 7*w_common + 20 - w_9 - w_10, we should minimize w_common.
    # We can fix w_9 and w_10 to their minimums to see the minimal w_common.
    # JR (a part of EJR) for voters 9 and 10 requires w_9 >= 2 and w_10 >= 2.
    w_9_ejr, w_10_ejr = 2, 2
    w_sum_N_ejr = committee_size_k - w_9_ejr - w_10_ejr
    
    print(f"To find the minimum satisfaction, we assume w_9={w_9_ejr} and w_10={w_10_ejr}.")
    print(f"This leaves w_common + w_pairs = {w_sum_N_ejr}.")

    # s(N,W) = 7*w_common + w_sum_N_ejr.
    # To satisfy EJR (s_j >= 8 for 5 voters), we need w_common + w_p_j >= 8.
    # With max w_p_j = 2, this means w_common + 2 >= 8, so w_common >= 6.
    w_common_min_ejr = 6
    print(f"To satisfy EJR with minimal cost, w_common must be at least {w_common_min_ejr}.")
    
    # Check if w_common=6 is feasible. w_pairs = 16 - 6 = 10.
    # To give 5 voters satisfaction >= 8, they need w_p_j >= 2.
    # This costs 5 * 2 = 10 places in w_pairs, which is exactly what we have. It is feasible.
    s_N_W2 = 7 * w_common_min_ejr + w_sum_N_ejr
    
    print(f"s(N, W_2) = 7 * {w_common_min_ejr} + {w_sum_N_ejr} = {s_N_W2}")
    print("-" * 30)

    # 4. Compute the final ratio
    print("--- Part 3: Final Ratio ---")
    ratio = s_N_W1 / s_N_W2
    
    print(f"The first satisfaction is s(N, W_1) = {s_N_W1}")
    print(f"The second satisfaction is s(N, W_2) = {s_N_W2}")
    print(f"The final equation is: {s_N_W1} / {s_N_W2}")
    
    # The final answer in decimal form
    print(f"The resulting ratio is {ratio}")
    
    # Final answer in the required format
    return ratio

if __name__ == "__main__":
    final_answer = solve_committee_election()
    print(f"<<<{final_answer}>>>")