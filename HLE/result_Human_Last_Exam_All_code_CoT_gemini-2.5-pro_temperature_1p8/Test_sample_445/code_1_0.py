def solve_probability_puzzle():
    """
    This function calculates the maximal probability p for Alice's guessing game
    based on the optimal strategy.
    """
    
    # The total number of boxes.
    n = 20
    
    # Alice's strategy involves selecting a random group of 'm' boxes.
    # The number of opened boxes 'k' is m-1.
    # The problem states 1 <= k <= 19, so 2 <= m <= 20.
    
    # The probability of success for a given group size 'm' is (m-2)/m.
    # This function is maximized when 'm' is as large as possible.
    optimal_m = 20
    
    # Calculate the numerator and denominator for the final probability.
    numerator = optimal_m - 2
    denominator = optimal_m
    
    # Calculate the final probability.
    probability = numerator / denominator
    
    print("The problem asks for the maximal guaranteed probability of success for Alice.")
    print("The optimal strategy is for Alice to form a group of m=20 boxes (all of them),")
    print("randomly choose 1 box to keep closed (the 'target'), and open the other 19.")
    print("Her probability of success P(m) follows the formula: P(m) = (m - 2) / m.")
    print(f"\nFor the optimal group size m = {optimal_m}:")
    
    # The final equation with each number explicitly shown.
    print(f"Maximal Probability = ({optimal_m} - 2) / {optimal_m}")
    print(f"Maximal Probability = {numerator} / {denominator}")
    print(f"Maximal Probability = {probability}")

# Run the function to display the result.
solve_probability_puzzle()