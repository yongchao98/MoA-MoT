def solve_hat_puzzle():
    """
    This function explains the solution to the hat puzzle and prints the results.
    """
    
    # The puzzle involves logical deduction, which can be reasoned through step-by-step.
    # The final state is determined by a process of elimination over three rounds of questions.
    
    # The reasoning leads to a specific type of arrangement for the 5 black and 4 white hats.
    # One such valid arrangement is BWBWBWBBW. In this arrangement, and all its rotations,
    # every white hat is situated between two black hats.
    
    # People with white hats are in a position to deduce their hat color in the third round.
    # People with black hats are not.
    
    num_black_hats = 5
    num_white_hats = 4
    
    # In the third round, everyone with a white hat can deduce their hat color.
    # The logic is that if they assumed their hat was black, it would create a scenario
    # that contradicts the fact that everyone answered "No" in Round 2.
    # This reasoning is available to all people with white hats simultaneously.
    
    num_people_saying_yes = num_white_hats
    
    # One possible distribution of hats. Note that any rotation or reflection of this
    # pattern is also a valid solution. For example, BBWBBWBWB is a reflection.
    hat_distribution = ['B', 'W', 'B', 'W', 'B', 'W', 'B', 'B', 'W']
    
    print("How many people replied 'Yes.'?")
    print(f"Number of people who replied 'Yes': {num_people_saying_yes}")
    print("\nHow are black and white hats distributed around the table?")
    print("One possible distribution is (read clockwise):")
    
    # We print the distribution with each hat color on a new line as part of the equation.
    # While it's not a mathematical equation, we can represent the setup this way.
    print("Person 1 hat = " + hat_distribution[0])
    print("Person 2 hat = " + hat_distribution[1])
    print("Person 3 hat = " + hat_distribution[2])
    print("Person 4 hat = " + hat_distribution[3])
    print("Person 5 hat = " + hat_distribution[4])
    print("Person 6 hat = " + hat_distribution[5])
    print("Person 7 hat = " + hat_distribution[6])
    print("Person 8 hat = " + hat_distribution[7])
    print("Person 9 hat = " + hat_distribution[8])

solve_hat_puzzle()
# The final answer for the number of people is 4.
# The final answer for the distribution is any pattern equivalent to BWBWBWBBW.