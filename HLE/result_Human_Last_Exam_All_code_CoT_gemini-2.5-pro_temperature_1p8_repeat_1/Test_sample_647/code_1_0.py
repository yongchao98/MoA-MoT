import math

def solve_committee_problem():
    """
    Solves the committee election problem by calculating s(N, W1) and s(N, W2)
    and then finding their ratio.
    """
    # Number of voters n and committee size k
    n = 10
    k = 20

    # Group of voters N = {1, ..., 8}
    N_size = 8

    # --- Part 1: Finding s(N, W_1) for a committee in the Core ---

    print("### Step 1: Analysis of the Core and finding s(N, W_1)\n")
    
    # In this profile, the approval scores are:
    # Candidates {1, ..., 8}: score 8 (each approved by 8 voters in N)
    # Candidates {9, ..., 32}: score 1 (each approved by one voter)
    
    # A committee in the "core" for this profile is an Approval-Voting-based committee.
    # It consists of the k=20 candidates with the highest approval scores.
    # This means any core committee W_core must contain all 8 candidates from {1..8}.
    # The remaining 20 - 8 = 12 members must be chosen from the other 24 candidates {9..32}.
    core_cands_from_1_8 = 8
    core_cands_from_9_32 = 12

    # The satisfaction of group N with a committee W is s(N,W) = sum_{i in N} |A(i) intersection W|.
    # For any core committee W_core, it contains {1..8}. So, for any voter i in N,
    # the intersection A(i) with W_core includes {1..8}.
    # s(N, W_core) = sum_{i=1 to 8} (|{1..8}| + |(A(i)\{1..8}) intersect W_core|)
    # s(N, W_core) = 8 * 8 + sum_{i=1 to 8} |P_i intersect W_12|, where P_i are the unique candidates
    # for voter i in N and W_12 is the set of 12 candidates from {9..32}.
    
    satisfaction_from_common_core = N_size * core_cands_from_1_8

    # The sets P_i for i in {1..8} are disjoint pairs: {9,10}, {11,12}, ..., {23,24}.
    # Their union, let's call it C_N_unique, is {9..24}, containing 16 candidates.
    # The other candidates are {25..32} (8 candidates), from voters 9 and 10.
    
    # To find W_1 that gives the LOWEST satisfaction to N, we must choose W_12
    # to minimize its intersection with C_N_unique.
    # We can choose all 8 candidates from {25..32}. These have zero overlap with C_N_unique.
    # We still need to choose 12 - 8 = 4 more candidates.
    # These 4 must come from C_N_unique.
    
    min_overlap_from_unique_cands = 4
    
    s_N_W1 = satisfaction_from_common_core + min_overlap_from_unique_cands

    print(f"A committee W_1 in the core with the lowest satisfaction for N includes:")
    print(f"- All {core_cands_from_1_8} common candidates {{1..8}}.")
    print(f"- 12 other candidates chosen to minimize overlap with N's other preferences.")
    print(f"The minimum additional satisfaction from these 12 candidates is {min_overlap_from_unique_cands}.")
    print(f"s(N, W_1) = {N_size} * {core_cands_from_1_8} + {min_overlap_from_unique_cands} = {s_N_W1}\n")


    # --- Part 2: Finding s(N, W_2) for an EJR committee ---

    print("### Step 2: Analysis of EJR and finding s(N, W_2)\n")

    # The EJR condition for this profile (n=10, k=20) implies:
    # 1. For any group S of size |S| >= l*n/k = l*0.5 that is l-cohesive,
    #    some voter i in S must have |A(i) intersect W| >= l.
    # This leads to two key constraints for any EJR committee W:
    # a) For every voter i in {1..10}, |A(i) intersect W| >= 2.
    #    (Derived with l=2, S={i}, |S|=1 >= 2*0.5=1).
    # b) For any 4 voters in N, at least one must be satisfied with >= 8 candidates.
    #    (Derived with l=8, S subset N, |S|=4 >= 8*0.5=4).
    #    This implies that the number of voters in N with satisfaction < 8 must be at most 3.
    #    So, at least 5 voters in N must have satisfaction >= 8.
    
    # Let W_2 be composed of:
    # - m candidates from C_N = {1..8}
    # - p candidates from C_N_unique = {9..24}
    # - q candidates from C_{9,10} = {25..32}
    # m+p+q=20.

    # To satisfy constraint (a) for voters 9 and 10, we need q >= 4. Let's pick q=4.
    # Then m+p = 16.
    
    # The satisfaction for N is s(N, W) = sum_{i=1..8} (m + |P_i intersect W_p|)
    # Since sum(|P_i intersect W_p|) = p, we have s(N, W) = 8*m + p.
    # Substituting p = 16-m, gives s(N,W) = 8*m + (16-m) = 7*m + 16.
    # To minimize s(N, W), we must minimize m.

    # We need to find the minimum m that allows satisfying EJR constraints.
    # At least 5 voters in N must have sat_i = m + |P_i intersect W_p| >= 8.
    # Since |P_i| = 2, the maximum |P_i intersect W_p| is 2.
    # So, we need m + 2 >= 8, which implies m >= 6.

    # Let's check if m=6 is feasible.
    # If m=6, then p=10.
    # We need to satisfy `6 + |P_i intersect W_p| >= 8` for 5 voters.
    # This means `|P_i intersect W_p| >= 2`. Since |P_i|=2, we must fully include their P_i sets.
    # For 5 voters, we need to include 5 sets of P_i, which requires p = 5 * 2 = 10.
    # This is feasible, as we have p=10 available slots.
    # For the other 3 voters in N, their satisfaction will be m + 0 = 6, which is >= 2.

    m_min = 6
    s_N_W2 = 7 * m_min + 16
    
    print("An EJR committee W_2 with the lowest satisfaction for N can be constructed with:")
    print(f"- m = {m_min} candidates from the common set {{1..8}}.")
    print(f"This is the minimum 'm' that can satisfy the EJR fairness constraints.")
    print(f"The total satisfaction for N is s(N, W_2) = 7*m + 16 = 7*{m_min} + 16 = {s_N_W2}\n")

    # --- Part 3: Final Ratio ---
    
    print("### Step 3: Final Calculation\n")
    
    ratio = s_N_W1 / s_N_W2
    
    print(f"The satisfaction for the core committee is s(N, W_1) = {s_N_W1}.")
    print(f"The satisfaction for the EJR committee is s(N, W_2) = {s_N_W2}.")
    print(f"The final equation is:")
    print(f"s(N, W_1) / s(N, W_2) = {s_N_W1} / {s_N_W2}")

    # Simplify fraction for clarity
    common_divisor = math.gcd(s_N_W1, s_N_W2)
    numerator = s_N_W1 // common_divisor
    denominator = s_N_W2 // common_divisor

    print(f"Simplified fraction: {numerator}/{denominator}")
    print(f"Decimal value: {ratio}")
    
    # Return the answer in the required format
    return ratio, f"<<<{numerator}/{denominator}>>>"

if __name__ == '__main__':
    result, formatted_answer = solve_committee_problem()
    # The final print to stdout is just the answer as requested.
    # The explanation is provided by the function's print statements.
    # For final submission format, we only need the last line.
    print(formatted_answer)
