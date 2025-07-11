import itertools

def solve_ejr_voter1_bounds():
    """
    Solves for the min and max number of candidates approved by Voter 1
    in a 5-member committee satisfying Extended Justified Representation (EJR).
    """

    # Step 1: Define Election Parameters
    ballots = {
        "Voter 1": {"x1", "x2", "x3", "y1", "z3"},
        "Voter 2": {"x1", "x2", "x3", "y2"},
        "Voter 3": {"x1", "x2", "x3", "y3"},
        "Voter 4": {"y4", "y5", "y6", "z1"},
        "Voter 5": {"y4", "y5", "y6", "z1"},
        "Voter 6": {"y4", "y5", "y6", "z2"},
        "Voter 7": {"y4", "y5", "y6", "z2"},
        "Voter 8": {"x1", "x2", "x3", "z1"},
        "Voter 9": {"x1", "x2", "x3", "z1"},
        "Voter 10": {"x1", "x2", "x3", "z1"},
    }
    num_voters = 10
    committee_size = 5
    quota = num_voters / committee_size

    print("Step 1: Election Parameters")
    print(f"Number of voters (n): {num_voters}")
    print(f"Committee size (k): {committee_size}")
    print(f"Quota (n/k): {int(quota)}\n")

    # Step 2: Identify Key Cohesive Groups and EJR Constraints
    print("Step 2: Applying EJR to identify constraints on the committee composition.")
    
    # Group A: Voters {V1, V2, V3, V8, V9, V10}
    group_A_voters = ["Voter 1", "Voter 2", "Voter 3", "Voter 8", "Voter 9", "Voter 10"]
    size_A = len(group_A_voters)
    intersection_A = set.intersection(*(ballots[v] for v in group_A_voters))
    cohesion_A = len(intersection_A)
    union_A = set.union(*(ballots[v] for v in group_A_voters))
    # EJR check: The group size 6 >= l * 2. With cohesion 3, we check for l=3.
    # 6 >= 3 * 2 is true. So this group must be 3-represented.
    constraint_A = 3
    print(f"Group A = {group_A_voters}")
    print(f"   - Size: {size_A}, Common Candidates: {intersection_A} (Cohesion: {cohesion_A})")
    print(f"   - EJR Condition (l=3): Size {size_A} >= 3 * quota ({int(quota)})")
    print(f"   - Constraint 1: Committee must contain at least {constraint_A} members from Group A's approved candidates: {union_A}\n")
    
    # Group B: Voters {V4, V5, V6, V7}
    group_B_voters = ["Voter 4", "Voter 5", "Voter 6", "Voter 7"]
    size_B = len(group_B_voters)
    intersection_B = set.intersection(*(ballots[v] for v in group_B_voters))
    cohesion_B = len(intersection_B)
    union_B = set.union(*(ballots[v] for v in group_B_voters))
    # EJR check: The group size 4 >= l * 2. With cohesion 3, we check for l=2.
    # 4 >= 2 * 2 is true. So this group must be 2-represented.
    constraint_B = 2
    print(f"Group B = {group_B_voters}")
    print(f"   - Size: {size_B}, Common Candidates: {intersection_B} (Cohesion: {cohesion_B})")
    print(f"   - EJR Condition (l=2): Size {size_B} >= 2 * quota ({int(quota)})")
    print(f"   - Constraint 2: Committee must contain at least {constraint_B} members from Group B's approved candidates: {union_B}\n")
    
    # Step 3: Characterize Committee Structure from Constraints
    print("Step 3: Characterizing committee structure based on constraints.")
    C_A_only = union_A - union_B
    C_B_only = union_B - union_A
    C_intersect = union_A & union_B
    
    # w_a = |W ∩ C_A_only|, w_b = |W ∩ C_B_only|, w_z1 = |W ∩ C_intersect|
    # w_a + w_z1 >= 3
    # w_b + w_z1 >= 2
    # w_a + w_b + w_z1 = 5
    print("This leads to 3 possible committee structures:")
    print("  1. Committee has z1, 2 from C_A_only, 2 from C_B_only.")
    print("  2. Committee has z1, 3 from C_A_only, 1 from C_B_only.")
    print("  3. Committee has no z1, 3 from C_A_only, 2 from C_B_only.\n")
    
    # Step 4: Calculate Min/Max Approved Candidates for Voter 1
    print("Step 4: Calculating min/max committee members approved by Voter 1 for each structure.")
    voter1_approved = ballots["Voter 1"]
    other_C_A_only = C_A_only - voter1_approved
    
    # Structure 1: w_a=2, w_b=2, w_z1=1
    s1_min = max(0, 2 - len(other_C_A_only))
    s1_max = min(2, len(voter1_approved))
    
    # Structure 2: w_a=3, w_b=1, w_z1=1
    s2_min = max(0, 3 - len(other_C_A_only))
    s2_max = min(3, len(voter1_approved))
    
    # Structure 3: w_a=3, w_b=2, w_z1=0
    s3_min = max(0, 3 - len(other_C_A_only))
    s3_max = min(3, len(voter1_approved))
    
    print(f"Voter 1's ballot: {voter1_approved}")
    print(f"Structure 1: Min approved = {s1_min}, Max approved = {s1_max}")
    print(f"Structure 2: Min approved = {s2_min}, Max approved = {s2_max}")
    print(f"Structure 3: Min approved = {s3_min}, Max approved = {s3_max}\n")

    # Step 5: Determine Overall Min/Max
    min_approved_v1 = min(s1_min, s2_min, s3_min)
    max_approved_v1 = max(s1_max, s2_max, s3_max)
    
    print("Step 5: Final Result")
    print(f"The overall minimum is min({s1_min}, {s2_min}, {s3_min}) = {min_approved_v1}")
    print(f"The overall maximum is max({s1_max}, {s2_max}, {s3_max}) = {max_approved_v1}\n")
    
    print(f"Conclusion: The minimum number of candidates in the committee that are approved by voter 1 is {min_approved_v1}.")
    print(f"Conclusion: The maximum number of candidates in the committee that are approved by voter 1 is {max_approved_v1}.")


solve_ejr_voter1_bounds()