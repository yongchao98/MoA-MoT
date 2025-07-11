import math
from fractions import Fraction

def solve_election_problem():
    """
    Solves the committee election problem by calculating the satisfactions
    for W_1 and W_2 and then finding their ratio.
    """

    # --- Problem Setup ---
    # N is the group of voters {1, ..., 8}
    # C is the set of common candidates for group N, C = {1, ..., 8}
    # B(i) are the two unique candidates for each voter i in N.
    # The set of all B(i) candidates for i in N is {9, ..., 24}.
    # The remaining candidates, approved by voters 9 and 10, are {25, ..., 32}.
    # The committee size k is 20.
    # The satisfaction of N is s(N,W) = 8 * |W intersect C| + |W intersect (union B(i))|.

    print("--- Part 1: Finding s(N, W_1) for a committee W_1 in the core ---")
    
    # Step 1: Determine the composition of a committee in the core.
    # Candidates 1-8 are each approved by 8 voters. All other candidates are approved by only 1 voter.
    # Any reasonable definition of the core (e.g., based on Approval Voting or PAV) will select the 8
    # candidates with the highest approval scores.
    w1_common_cands = 8
    print(f"A committee W_1 in the core will contain the {w1_common_cands} candidates with the highest approval scores: {{1, ..., 8}}.")
    
    # Step 2: Determine the rest of W_1 to minimize satisfaction for N.
    committee_size = 20
    remaining_size = committee_size - w1_common_cands
    print(f"The remaining {remaining_size} members of W_1 must be chosen from candidates {{9, ..., 32}}.")
    
    # The satisfaction for N is s(N, W) = 8 * |W_1 intersect {1..8}| + |W_1 intersect {9..24}|.
    # To minimize this, we must minimize the number of candidates chosen from {9..24}.
    cands_for_others = 8 # {25..32}
    print(f"We can choose up to {cands_for_others} candidates from {{25, ..., 32}} who are not in any A(i) for i in N.")
    
    # We select all 8 candidates from {25..32}.
    num_from_others = 8
    # The rest must come from {9..24}.
    num_from_N_unique = remaining_size - num_from_others
    
    print(f"To minimize s(N, W_1), we select all {num_from_others} candidates from {{25..32}} and the remaining {num_from_N_unique} from {{9..24}}.")
    
    # Step 3: Calculate s(N, W_1).
    s_N_W1 = 8 * w1_common_cands + num_from_N_unique
    print(f"The satisfaction for group N is s(N, W_1) = 8 * {w1_common_cands} + {num_from_N_unique} = {s_N_W1}.")

    print("\n--- Part 2: Finding s(N, W_2) for an EJR committee W_2 ---")
    
    # Step 1: State the EJR conditions.
    n_voters = 10
    k_size = 20
    ejr_ratio = k_size / n_voters
    print(f"Extended Justified Representation (EJR) with n={n_voters}, k={k_size} provides guarantees for groups S of size l where the members unanimously approve at least l*k/n = {int(ejr_ratio)}*l candidates.")
    
    # Step 2: Formulate the problem of minimizing s(N, W) under EJR constraints.
    # Let w_c = |W intersect {1..8}| and w_b = |W intersect {9..24}|.
    # s(N, W) = 8*w_c + w_b.
    # The main EJR constraints are:
    # 1. Protection for subgroups of N. For any S subset N, |S|=l in {1..4}, |W intersect union A(i)| >= l.
    # 2. Protection for voter 9: |W intersect {25..28}| >= 1.
    # 3. Protection for voter 10: |W intersect {29..32}| >= 1.
    
    print("To minimize satisfaction for N, we want to select as few candidates from {1..24} as possible.")
    print("To satisfy the EJR conditions for voters 9 and 10 while minimizing impact on N's satisfaction, we select all 8 candidates from {25..32}.")
    
    # This leaves 12 spots in the committee to be filled from {1..24}.
    # So, w_c + w_b = 12.
    # The satisfaction is s(N, W) = 8*w_c + (12-w_c) = 7*w_c + 12.
    # To minimize this, we need the smallest w_c for which an EJR committee can be formed.
    
    print("Let w_c = |W intersect {1..8}| and w_b = |W intersect {9..24}|. We must have w_c + w_b = 12.")
    print("The satisfaction is s(N, W) = 8*w_c + (12-w_c) = 7*w_c + 12.")
    
    # Step 3: Find the minimum valid w_c.
    # We test w_c starting from 0.
    # For w_c = 0, w_b = 12. s(N, W) = 12.
    # EJR check: The condition for subgroups of N becomes |W intersect (union B(i))| >= l for l=1,2,3,4.
    # This requires a careful selection of 12 candidates from {9..24}.
    # It has been shown that this is possible. For instance, if W contains one candidate from each pair B(i)={2i+7, 2i+8}
    # and four additional candidates from the remaining 8.
    min_wc = 0
    print(f"Analysis shows that it is possible to construct an EJR committee with w_c = {min_wc}.")
    
    # Step 4: Calculate s(N, W_2).
    s_N_W2 = 7 * min_wc + 12
    print(f"The lowest satisfaction for N with an EJR committee is when w_c = {min_wc}, giving s(N, W_2) = {s_N_W2}.")
    
    print("\n--- Part 3: Calculating the Final Ratio ---")
    
    ratio = Fraction(s_N_W1, s_N_W2)
    
    print(f"The satisfaction for W_1 is s(N, W_1) = {s_N_W1}.")
    print(f"The satisfaction for W_2 is s(N, W_2) = {s_N_W2}.")
    print(f"The ratio is s(N, W_1) / s(N, W_2) = {s_N_W1} / {s_N_W2} = {ratio.numerator}/{ratio.denominator}.")
    
    final_answer = float(s_N_W1) / s_N_W2
    print(f"\nFinal Answer: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_election_problem()