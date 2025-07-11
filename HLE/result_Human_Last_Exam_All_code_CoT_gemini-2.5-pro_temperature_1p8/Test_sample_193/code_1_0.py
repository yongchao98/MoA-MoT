import math

def calculate_optimal_poker_strategy():
    """
    This function calculates the optimal poker strategy for Hero with AA based on GTO principles
    for the given scenario.
    """

    # Define game parameters
    pot = 10
    stack = 1000

    # The problem boils down to finding the bet size 'B' that maximizes Hero's profit with AA.
    # In a GTO equilibrium for this scenario, the profit for Hero when value-betting size 'B' is
    # derived from the conditions that make both Hero indifferent to bluffing and Villain indifferent to calling.
    # The resulting profit function for a value bet of size B is:
    # Profit(B) = Pot + (Pot * B) / (Pot + B)
    #
    # We can analyze this function to find the optimal bet size.
    # The derivative of this function with respect to B is always positive, which means
    # the profit increases as the bet size B increases.
    # Therefore, to maximize profit, Hero should choose the largest possible bet size.
    optimal_bet_size = stack

    # The profit from betting this size is greater than the profit from checking (which is just the pot).
    # Profit(Bet 1000) = 10 + (10 * 1000) / (10 + 1000) = 10 + 9.9 = 19.9
    # Profit(Check) = 10
    # Since betting is more profitable, Hero should always bet with AA and never check.
    # The frequency of this action is 100%.
    strategy_frequency = 100

    # Per the user's request, we round the sizing and percentage to the nearest even number.
    def round_to_even(n):
        return round(n / 2.0) * 2

    final_sizing = round_to_even(optimal_bet_size)
    final_frequency = round_to_even(strategy_frequency)

    # Output the final strategy string.
    # The action is to BET, with the calculated sizing and 100% frequency.
    action = "BET"
    
    # We print the final equation which consists of the Action, Sizing, and Frequency.
    print(f"{action} {final_sizing} {final_frequency}%")

calculate_optimal_poker_strategy()