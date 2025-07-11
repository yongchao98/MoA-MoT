def solve_puzzle():
    """
    Calculates the maximal probability p for Alice's guessing game.
    """
    # Total number of boxes in the sequence.
    total_boxes = 20

    # Alice's strategy is to open 19 boxes and guess the content of the last one.
    # The guess is that the number is between the minimum and maximum of the opened boxes.
    # This strategy fails only if the box left closed contains an "extreme" value.
    # The extreme values are the global minimum and the global maximum of all 20 numbers.
    losing_cases = 2

    # The number of winning cases is the total number of boxes minus the number of losing cases.
    winning_cases = total_boxes - losing_cases

    # The probability of success is the ratio of winning cases to the total number of boxes
    # that Alice could choose to leave closed.
    probability = winning_cases / total_boxes

    print("The problem asks for the maximal guaranteed probability of success for Alice.")
    print("The optimal strategy yielding this probability is as follows:")
    print("1. Alice randomly selects 1 box to leave closed.")
    print(f"2. She opens the other {total_boxes - 1} boxes.")
    print("3. She guesses the number in the closed box is between the minimum and maximum of the numbers she observed.")
    print("\nThis strategy succeeds if the closed box does not contain the global minimum or global maximum value.")
    
    print("\nLet's calculate the probability:")
    print(f"Total number of boxes Alice can choose to leave closed: {total_boxes}")
    print(f"Number of 'losing' choices (the box with the min or max value): {losing_cases}")
    print(f"Number of 'winning' choices is the total minus the losing choices.")
    print(f"Winning choices = {total_boxes} - {losing_cases} = {winning_cases}")
    
    print("\nThe probability is the ratio of winning choices to total choices:")
    print(f"p = {winning_cases} / {total_boxes}")
    print(f"p = {probability}")
    print("\nIn fraction form, the maximal probability p is 9/10.")

solve_puzzle()