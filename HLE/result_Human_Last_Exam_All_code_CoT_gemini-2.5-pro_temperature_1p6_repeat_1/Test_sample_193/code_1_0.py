import math

def solve_poker_strategy():
    """
    This function calculates and prints the optimal strategy for AA in the given scenario.
    The logic is derived in the text above, concluding that the optimal strategy for AA
    is to bet the maximum possible amount, which is all-in.
    """
    # Player's action: Always bet with the nut hand (AA).
    action = "BET"

    # Bet Sizing: The EV-maximizing strategy is to bet the largest possible amount.
    # Stacks are $1000 behind.
    bet_sizing = 1000

    # Frequency: Since betting is always superior to checking for AA in this spot,
    # it's a pure strategy, played 100% of the time.
    frequency_pct = 100

    # The problem asks to round to the nearest even number.
    # Sizing 1000 is even.
    # Percentage 100 is even.
    final_sizing = int(bet_sizing)
    final_frequency = int(frequency_pct)
    
    # Print the final answer in the specified format.
    print(f"{action} {final_sizing} {final_frequency}%")

solve_poker_strategy()