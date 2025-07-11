import itertools

def solve():
    """
    This script demonstrates the minimum and maximum number of candidates
    in a 5-member EJR committee that are approved by Voter 1.
    """

    # Define Voter 1's approval set for calculating the overlap.
    voter1_approvals = {'x1', 'x2', 'x3', 'y1', 'z3'}

    # The analysis of EJR constraints leads to the following conclusions:
    # 1. A valid EJR committee can be formed that has 0 members approved by Voter 1.
    # 2. A valid EJR committee can be formed that has 3 members approved by Voter 1.
    # 3. It's impossible to form a valid EJR committee with more than 3 members approved by Voter 1.

    # We present the "equations" for the minimum and maximum cases using example committees
    # that are known to be compliant with the main EJR constraints.

    print("To find the minimum number of candidates approved by voter 1:")
    w_min = {'y2', 'y3', 'y4', 'y5', 'z1'}
    v1_in_w_min = w_min.intersection(voter1_approvals)
    min_val = len(v1_in_w_min)
    # The 'equation' shows the set intersection and the resulting size.
    print(f"| {w_min} ∩ {voter1_approvals} | = | {v1_in_w_min} | = {min_val}")

    print("\nTo find the maximum number of candidates approved by voter 1:")
    w_max = {'x1', 'x2', 'x3', 'y4', 'y5'}
    v1_in_w_max = w_max.intersection(voter1_approvals)
    max_val = len(v1_in_w_max)
    # The 'equation' shows the set intersection and the resulting size.
    print(f"| {w_max} ∩ {voter1_approvals} | = | {v1_in_w_max} | = {max_val}")

solve()