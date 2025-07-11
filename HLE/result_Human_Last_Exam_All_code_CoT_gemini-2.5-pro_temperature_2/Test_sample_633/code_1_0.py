import collections
from itertools import combinations

def solve_ejr_satisfaction():
    """
    This script analyzes a voter profile to determine the min and max
    satisfaction for Voter 1 in any committee satisfying Extended
    Justified Representation (EJR).
    """
    
    # --- Setup: Define voter ballots, committee size, and quota ---
    ballots = {
        'Voter 1': {'x1', 'x2', 'x3', 'y1', 'z3'},
        'Voter 2': {'x1', 'x2', 'x3', 'y2'},
        'Voter 3': {'x1', 'x2', 'x3', 'y3'},
        'Voter 4': {'y4', 'y5', 'y6', 'z1'},
        'Voter 5': {'y4', 'y5', 'y6', 'z1'},
        'Voter 6': {'y4', 'y5', 'y6', 'z2'},
        'Voter 7': {'y4', 'y5', 'y6', 'z2'},
        'Voter 8': {'x1', 'x2', 'x3', 'z1'},
        'Voter 9': {'x1', 'x2', 'x3', 'z1'},
        'Voter 10': {'x1', 'x2', 'x3', 'z1'},
    }
    committee_size_k = 5
    num_voters_n = len(ballots)
    quota = num_voters_n / committee_size_k

    print(f"Number of voters (n): {num_voters_n}")
    print(f"Committee size (k): {committee_size_k}")
    print(f"Hare quota (n/k): {num_voters_n} / {committee_size_k} = {quota}")
    print("-" * 50)

    # --- Step 1: Analyze cohesive groups and apply EJR rules ---
    print("Step 1: Analyzing cohesive voter groups for EJR constraints...\n")
    
    # Group A: Voters unanimously approving {x1, x2, x3}
    group_A_candidates = {'x1', 'x2', 'x3'}
    group_A_voters = [v for v, b in ballots.items() if group_A_candidates.issubset(b)]
    j_A = len(group_A_voters) / quota
    print(f"A group of {len(group_A_voters)} voters all approve {len(group_A_candidates)} candidates ({group_A_candidates}).")
    print(f"Their size relative to the quota is j = {len(group_A_voters)} / {quota} = {int(j_A)}.")
    print(f"EJR requires the committee to contain at least {int(j_A)} members from this set.")
    print(f"=> This implies all of {list(group_A_candidates)} must be in any EJR committee.\n")
    
    # Group B: Voters unanimously approving {y4, y5, y6}
    group_B_candidates = {'y4', 'y5', 'y6'}
    group_B_voters = [v for v, b in ballots.items() if group_B_candidates.issubset(b)]
    j_B = len(group_B_voters) / quota
    print(f"A group of {len(group_B_voters)} voters all approve {len(group_B_candidates)} candidates ({group_B_candidates}).")
    print(f"Their size relative to the quota is j = {len(group_B_voters)} / {quota} = {int(j_B)}.")
    print(f"EJR requires the committee to contain at least {int(j_B)} members from this set.\n")
    print("-" * 50)

    # --- Step 2: Characterize and construct all possible EJR committees ---
    print("Step 2: Constructing all possible EJR committees based on constraints...\n")
    
    must_include = {'x1', 'x2', 'x3'}
    must_include_from_pool = {'y4', 'y5', 'y6'}
    must_include_from_count = 2 # from Group B's constraint
    
    remaining_slots = committee_size_k - len(must_include)
    print(f"A committee must include {must_include}, leaving {remaining_slots} slots to fill.")
    print(f"These {remaining_slots} slots must be filled by {must_include_from_count} candidates from {must_include_from_pool}.")
    
    possible_committees = []
    additional_choices = combinations(must_include_from_pool, must_include_from_count)
    for choice in additional_choices:
        committee = must_include.union(set(choice))
        possible_committees.append(frozenset(committee))
        
    print(f"\nThis results in {len(possible_committees)} possible EJR committees:")
    for c in possible_committees:
        print(f"  - {sorted(list(c))}")
    print("-" * 50)
    
    # --- Step 3: Calculate satisfaction for Voter 1 for each committee ---
    print("Step 3: Calculating satisfaction for Voter 1 for each possible committee...\n")
    voter_1_ballot = ballots['Voter 1']
    print(f"Voter 1's approval ballot is: {sorted(list(voter_1_ballot))}\n")

    satisfaction_scores = []
    for committee in possible_committees:
        approved_in_committee = committee.intersection(voter_1_ballot)
        score = len(approved_in_committee)
        satisfaction_scores.append(score)
        print(f"For committee {sorted(list(committee))}:")
        print(f"  Approved members: {sorted(list(approved_in_committee))}")
        # The line below shows the numbers in the "final equation" for this committee
        print(f"  Satisfaction = |Voter 1 ballot ∩ Committee| = {score}\n")

    # --- Step 4: Determine the minimum and maximum satisfaction ---
    min_satisfaction = min(satisfaction_scores) if satisfaction_scores else "N/A"
    max_satisfaction = max(satisfaction_scores) if satisfaction_scores else "N/A"

    print("-" * 50)
    print(f"Final Answer:")
    print(f"The minimum number of candidates approved by Voter 1 is: {min_satisfaction}")
    print(f"The maximum number of candidates approved by Voter 1 is: {max_satisfaction}")

solve_ejr_satisfaction()