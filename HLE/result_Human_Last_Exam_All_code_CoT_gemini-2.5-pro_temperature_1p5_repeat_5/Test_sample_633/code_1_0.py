import itertools

def solve_ejr_bounds():
    """
    This script finds the minimum and maximum number of candidates
    approved by Voter 1 in a 5-member committee that satisfies
    Extended Justified Representation (EJR).
    """
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

    n = len(voters)
    k = 5
    q = n // k + 1

    all_candidates = set.union(*voters.values())
    voter_list = list(voters.values())
    voter1_approved_set = voter_list[0]

    memo_group_props = {}

    def get_group_properties(group_indices):
        """Calculates and caches properties of a voter group."""
        group_indices_tuple = tuple(sorted(group_indices))
        if group_indices_tuple in memo_group_props:
            return memo_group_props[group_indices_tuple]

        group_voter_sets = [voter_list[i] for i in group_indices]
        
        intersection = group_voter_sets[0].copy()
        for i in range(1, len(group_voter_sets)):
            intersection.intersection_update(group_voter_sets[i])
            
        union = set.union(*group_voter_sets)
        
        memo_group_props[group_indices_tuple] = (intersection, union)
        return intersection, union

    def is_ejr(committee):
        """Checks if a committee satisfies Extended Justified Representation."""
        voter_indices = range(n)
        for l in range(1, k + 1):
            min_group_size = l * q
            if min_group_size > n:
                break
            
            for r in range(min_group_size, n + 1):
                for group_indices in itertools.combinations(voter_indices, r):
                    intersection, union = get_group_properties(group_indices)
                    
                    if len(intersection) >= l:  # The group is l-cohesive
                        # Check EJR condition
                        if len(committee.intersection(union)) < l:
                            return False
        return True

    min_approved = k + 1
    max_approved = -1
    
    print("Starting analysis...")
    print(f"Number of voters (n) = {n}")
    print(f"Committee size (k) = {k}")
    print(f"Droop quota (q) = floor(n/k) + 1 = floor({n}/{k}) + 1 = {q}")
    print(f"Total number of unique candidates: {len(all_candidates)}")
    print(f"Total possible committees of size {k}: {len(list(itertools.combinations(all_candidates, k)))}")
    print("-" * 20)

    # Iterate through all possible committees
    for committee_tuple in itertools.combinations(all_candidates, k):
        committee = set(committee_tuple)
        if is_ejr(committee):
            approved_count = len(committee.intersection(voter1_approved_set))
            min_approved = min(min_approved, approved_count)
            max_approved = max(max_approved, approved_count)

    print("Analysis Complete.")
    print("\n--- Final Result ---")
    if max_approved == -1:
        print("No EJR committee was found.")
    else:
        print(f"The minimum number of candidates in an EJR committee approved by voter 1 is: {min_approved}")
        print(f"The maximum number of candidates in an EJR committee approved by voter 1 is: {max_approved}")
        print("\nThe final equation for the minimum number of approved candidates for Voter 1 is:")
        print(f"min |W ∩ A_1| = {min_approved}")
        print("\nThe final equation for the maximum number of approved candidates for Voter 1 is:")
        print(f"max |W ∩ A_1| = {max_approved}")

solve_ejr_bounds()