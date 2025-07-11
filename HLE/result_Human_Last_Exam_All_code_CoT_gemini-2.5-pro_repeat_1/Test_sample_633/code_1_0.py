import itertools

def solve():
    """
    This script finds the minimum and maximum number of candidates approved by Voter 1
    in a 5-member committee that satisfies Extended Justified Representation (EJR).
    """
    
    # Define the preference profile of the 10 voters
    voters_str = {
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
    
    voters_list = list(voters_str.values())
    n = len(voters_list)  # Number of voters
    k = 5  # Committee size

    # Collect all unique candidates
    all_candidates = sorted(list(set.union(*voters_list)))

    def is_ejr(committee, ballots, n):
        """
        Checks if a given committee satisfies Extended Justified Representation (EJR).
        A committee is EJR if for any l-cohesive group of l voters, at least one
        voter in the group has at least l of their approved candidates in the committee.
        """
        committee_set = set(committee)
        voter_indices = range(n)
        
        for l in range(1, n + 1):
            # Iterate through all subgroups of voters of size l
            for group_indices in itertools.combinations(voter_indices, l):
                group_ballots = [ballots[i] for i in group_indices]
                
                # Find the set of commonly approved candidates for the group
                intersection = group_ballots[0].copy()
                for i in range(1, l):
                    intersection.intersection_update(group_ballots[i])
                
                # Check if the group is l-cohesive
                if len(intersection) >= l:
                    # If cohesive, check if the EJR condition holds
                    is_satisfied = False
                    for idx in group_indices:
                        voter_ballot = ballots[idx]
                        satisfaction = len(committee_set.intersection(voter_ballot))
                        if satisfaction >= l:
                            is_satisfied = True
                            break
                    
                    if not is_satisfied:
                        # Found a violation of EJR
                        return False
        return True

    # Initialize min and max approved counts
    min_approved = float('inf')
    max_approved = float('-inf')
    found_ejr_committee = False
    voter1_ballot = voters_str['Voter 1']

    # Iterate through all possible committees of size k
    for committee in itertools.combinations(all_candidates, k):
        if is_ejr(committee, voters_list, n):
            found_ejr_committee = True
            
            # Calculate how many candidates Voter 1 approves in this committee
            approved_by_v1 = len(set(committee).intersection(voter1_ballot))
            
            min_approved = min(min_approved, approved_by_v1)
            max_approved = max(max_approved, approved_by_v1)
            
    if found_ejr_committee:
        print(f"For an EJR committee of size {k}:")
        print(f"The minimum number of candidates approved by voter 1 is: {min_approved}")
        print(f"The maximum number of candidates approved by voter 1 is: {max_approved}")
    else:
        print(f"No EJR committee of size {k} was found.")

solve()