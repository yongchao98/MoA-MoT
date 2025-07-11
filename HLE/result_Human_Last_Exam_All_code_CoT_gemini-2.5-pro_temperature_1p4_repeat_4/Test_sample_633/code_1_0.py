import itertools

def solve():
    """
    Solves the EJR committee problem by finding the min and max overlap for Voter 1.
    """
    voter_ballots = [
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

    k = 5
    n = 10
    
    # All unique candidates in the election
    all_candidates = set.union(*voter_ballots)
    
    # Voter 1's ballot
    voter_1_ballot = voter_ballots[0]
    
    # From our analysis, we identified two key cohesive groups and their unions
    # Group X: V{1,2,3,8,9,10}, l=3
    group_x_indices = [0, 1, 2, 7, 8, 9]
    union_x = set.union(*(voter_ballots[i] for i in group_x_indices))
    l_x = 3
    
    # Group Y: V{4,5,6,7}, l=2
    group_y_indices = [3, 4, 5, 6]
    union_y = set.union(*(voter_ballots[i] for i in group_y_indices))
    l_y = 2
    
    min_approved = float('inf')
    max_approved = float('-inf')
    ejr_committee_found = False

    # Generate all possible committees of size k
    for committee_tuple in itertools.combinations(all_candidates, k):
        committee = set(committee_tuple)
        
        # Check if the committee satisfies the EJR constraints
        # Constraint 1 for Group X
        if len(committee.intersection(union_x)) < l_x:
            continue
            
        # Constraint 2 for Group Y
        if len(committee.intersection(union_y)) < l_y:
            continue

        # If we reach here, the committee is a valid EJR committee (based on our key constraints)
        ejr_committee_found = True
        
        # Calculate the number of approved candidates for Voter 1
        approved_count = len(committee.intersection(voter_1_ballot))
        
        if approved_count < min_approved:
            min_approved = approved_count
            
        if approved_count > max_approved:
            max_approved = approved_count

    if ejr_committee_found:
        print(f"The number of candidates in the committee is {k}.")
        print(f"The minimum number of candidates approved by voter 1 is {min_approved}.")
        print(f"The maximum number of candidates approved by voter 1 is {max_approved}.")
    else:
        print("No EJR committee was found based on the checked constraints.")

solve()