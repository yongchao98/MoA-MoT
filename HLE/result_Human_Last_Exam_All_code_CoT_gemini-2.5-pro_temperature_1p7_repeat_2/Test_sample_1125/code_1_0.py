def solve_puzzle():
    """
    This function determines and explains the solution to the hat puzzle.
    """
    total_members = 12
    # The strategy involves sacrificing two members to act as anchors or reference points.
    # The remaining members can use these anchors to determine their own numbers.
    num_anchors = 2
    
    # The number of people who are guaranteed to determine their number is the total
    # number of members minus the number of anchors needed for the strategy to work.
    guaranteed_winners = total_members - num_anchors
    
    print(f"The total number of team members is {total_members}.")
    print(f"A number of {num_anchors} members are designated as 'anchors' to create a frame of reference.")
    print("The remaining members are the 'solvers'.")
    print(f"Number of solvers = Total Members - Number of Anchors")
    print(f"Number of solvers = {total_members} - {num_anchors}")
    print(f"Number of solvers = {guaranteed_winners}")
    print("\nEach of the 10 solvers can determine their hat number with certainty.")
    print("The 2 anchors cannot, as all information is symmetrical with respect to them.")
    print("Therefore, the largest possible value of N is 10.")
    

solve_puzzle()