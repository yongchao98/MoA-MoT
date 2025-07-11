def solve_committee_problem():
    """
    Solves the EJR committee problem based on logical deduction.
    """
    # 1. Define the voter profile
    voter_profiles = {
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
    
    # 2. Define problem parameters
    num_voters = len(voter_profiles)
    committee_size = 5
    voter_share = num_voters / committee_size
    
    print(f"Number of voters (n): {num_voters}")
    print(f"Committee size (k): {committee_size}")
    print(f"Required voters per seat (n/k): {voter_share}")
    print("-" * 30)

    # 3. Analyze key voter groups and apply EJR principles
    
    # Group X analysis
    group_x_voters = [1, 2, 3, 8, 9, 10]
    group_x_size = len(group_x_voters)
    group_x_entitlement = group_x_size / voter_share
    group_x_common_candidates = set.intersection(*[voter_profiles[f'Voter {i}'] for i in group_x_voters])
    
    print("Analysis of Group X (Voters {1, 2, 3, 8, 9, 10}):")
    print(f"  - Size: {group_x_size}")
    print(f"  - Entitlement: {group_x_size} / {voter_share} = {int(group_x_entitlement)} seats")
    print(f"  - Common candidates: {group_x_common_candidates}")
    print("  - EJR Conclusion: This group is perfectly cohesive on 3 candidates and is entitled to 3 seats.")
    print("    Therefore, these 3 candidates must be in any EJR committee.")
    
    mandatory_members = group_x_common_candidates
    print(f"Mandatory Committee Members: {sorted(list(mandatory_members))}")
    print("-" * 30)

    # Group Y analysis
    group_y_voters = [4, 5, 6, 7]
    group_y_size = len(group_y_voters)
    group_y_entitlement = group_y_size / voter_share
    group_y_approved_pool = set.union(*[voter_profiles[f'Voter {i}'] for i in group_y_voters])

    print("Analysis of Group Y (Voters {4, 5, 6, 7}):")
    print(f"  - Size: {group_y_size}")
    print(f"  - Entitlement: {group_y_size} / {voter_share} = {int(group_y_entitlement)} seats")
    print("  - EJR Conclusion: This group is entitled to 2 seats, which must come from their approved candidates.")

    remaining_seats_pool = group_y_approved_pool
    print(f"Pool for remaining 2 seats: {sorted(list(remaining_seats_pool))}")
    print("-" * 30)

    # 4. Determine approvals for Voter 1
    voter1_ballot = voter_profiles['Voter 1']
    print(f"Calculating approvals for Voter 1:")
    print(f"Voter 1's ballot: {sorted(list(voter1_ballot))}")
    
    # Intersection with mandatory members
    approved_mandatory = voter1_ballot.intersection(mandatory_members)
    num_approved_mandatory = len(approved_mandatory)

    # Intersection with the pool for remaining members
    approved_from_pool = voter1_ballot.intersection(remaining_seats_pool)
    num_approved_from_pool = len(approved_from_pool)

    print(f"Voter 1's candidates from the mandatory set: {sorted(list(approved_mandatory))}")
    print(f"Voter 1's candidates from the remaining pool: {sorted(list(approved_from_pool))}")
    
    total_approved = num_approved_mandatory + num_approved_from_pool
    
    print("\nFinal Equation:")
    print(f"| Voter_1_Approvals | = | A_1 \cap W |")
    print(f"                    = | A_1 \cap (Mandatory \cup Remaining) |")
    print(f"                    = | (A_1 \cap Mandatory) \cup (A_1 \cap Remaining) |")
    # Print with actual sets
    v1_ballot_str = str(sorted(list(voter1_ballot)))
    mandatory_str = str(sorted(list(mandatory_members)))
    remaining_pool_str = str(sorted(list(remaining_seats_pool)))
    print(f"                    = | ({v1_ballot_str} \cap {mandatory_str}) \cup ({v1_ballot_str} \cap {remaining_pool_str}) |")
    print(f"                    = | {str(sorted(list(approved_mandatory)))} \cup {str(sorted(list(approved_from_pool)))} |")
    
    print(f"                    = {num_approved_mandatory} + {num_approved_from_pool}")
    print(f"                    = {total_approved}")
    
    print("-" * 30)
    print(f"Conclusion: The number of candidates approved by Voter 1 in any EJR committee is exactly {total_approved}.")
    print(f"Minimum number of approved candidates: {total_approved}")
    print(f"Maximum number of approved candidates: {total_approved}")

solve_committee_problem()
<<<E>>>