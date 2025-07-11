import itertools

def solve():
    """
    Solves the EJR committee problem by systematically checking all possible committees.
    """
    # 1. Define the Problem Data
    voters = {
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
    k = 5  # committee size
    n = len(voters) # number of voters
    voter_list = list(voters.items())

    all_candidates = set()
    for approved_set in voters.values():
        all_candidates.update(approved_set)
    
    # 2. Implement the EJR checker function
    def is_ejr(committee, voters_list, n, k):
        """
        Checks if a committee satisfies Extended Justified Representation (EJR).
        """
        num_voters = len(voters_list)
        # Iterate through all possible subgroup sizes l
        for l in range(1, num_voters + 1):
            # Iterate through all subgroups of voters of size l
            for voter_subgroup_indices in itertools.combinations(range(num_voters), l):
                
                subgroup_approvals = [voters_list[i][1] for i in voter_subgroup_indices]
                
                # Check for cohesiveness
                if not subgroup_approvals:
                    continue
                intersection = set.intersection(*subgroup_approvals)
                
                if not intersection:  # If not cohesive, this group poses no EJR constraint
                    continue

                # If cohesive, check the EJR condition
                total_representation = 0
                for i in voter_subgroup_indices:
                    voter_approval_set = voters_list[i][1]
                    representation = len(voter_approval_set.intersection(committee))
                    total_representation += representation
                    
                ejr_threshold = (l * k) // n # floor(l * k / n)
                
                if total_representation < ejr_threshold:
                    return False # EJR condition failed
        
        return True # EJR holds for all cohesive subgroups

    # 3. & 4. & 5. Iterate through committees, check EJR, and find min/max
    min_approved = float('inf')
    max_approved = float('-inf')
    min_committee_example = None
    max_committee_example = None
    ejr_committees_found = 0
    total_committees = 0

    for committee_tuple in itertools.combinations(all_candidates, k):
        total_committees += 1
        committee = set(committee_tuple)
        
        if is_ejr(committee, voter_list, n, k):
            ejr_committees_found += 1
            voter1_approvals = voters['Voter 1']
            num_approved_by_voter1 = len(voter1_approvals.intersection(committee))
            
            if num_approved_by_voter1 < min_approved:
                min_approved = num_approved_by_voter1
                min_committee_example = committee
            
            if num_approved_by_voter1 > max_approved:
                max_approved = num_approved_by_voter1
                max_committee_example = committee

    # 6. Print the final results
    print(f"Searching through {total_committees} possible committees for a valid EJR committee...")
    if ejr_committees_found > 0:
        print(f"Found {ejr_committees_found} EJR-compliant committees.\n")
        
        print(f"--- Final Answer ---")
        print(f"The minimum number of candidates approved by Voter 1 is: {min_approved}")
        print(f"The maximum number of candidates approved by Voter 1 is: {max_approved}\n")
        
        print("--- Example for Minimum ---")
        v1_rep_min = len(voters['Voter 1'].intersection(min_committee_example))
        print(f"An example committee that yields the minimum is W = {sorted(list(min_committee_example))}")
        print(f"Voter 1's representation = |{sorted(list(voters['Voter 1']))} ∩ {sorted(list(min_committee_example))}| = {v1_rep_min}")

        print("\n--- Example for Maximum ---")
        v1_rep_max = len(voters['Voter 1'].intersection(max_committee_example))
        print(f"An example committee that yields the maximum is W = {sorted(list(max_committee_example))}")
        print(f"Voter 1's representation = |{sorted(list(voters['Voter 1']))} ∩ {sorted(list(max_committee_example))}| = {v1_rep_max}")
    else:
        print("No EJR-compliant committees were found.")

solve()