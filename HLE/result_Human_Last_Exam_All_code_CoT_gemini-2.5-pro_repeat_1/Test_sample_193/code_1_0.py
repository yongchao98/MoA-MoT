import math

def round_to_nearest_even(n):
    """Rounds a number to the nearest even integer."""
    return int(round(n / 2) * 2)

def solve_poker_strategy():
    """
    Solves for the optimal poker strategy in the given scenario.
    """
    # Game parameters
    pot = 10
    stack = 1000

    # Step 1: Determine the optimal bet size.
    # To maximize the value of our bluffs, we need to bet as large as possible.
    # The optimal bet size is therefore all-in.
    bet_sizing = stack

    # Step 2: Determine the strategy for AA.
    # In the GTO equilibrium for this game, the hand AA is made indifferent
    # between betting and checking. This means any mixing frequency is valid.
    # A robust and canonical strategy is to play a 50/50 mix.
    bet_percentage = 50
    check_percentage = 50

    # Step 3: Round the results as per the user's request.
    final_sizing = round_to_nearest_even(bet_sizing)
    final_bet_percentage = round_to_nearest_even(bet_percentage)
    final_check_percentage = round_to_nearest_even(check_percentage)

    # Step 4: Format the output string, sorting actions by money put in.
    # The actions are BET and CHECK.
    bet_action = f"BET {final_sizing} {final_bet_percentage}%"
    check_action = f"CHECK {final_check_percentage}%"

    # The final strategy is a mix of betting and checking.
    final_strategy = f"{bet_action}, {check_action}"

    print(final_strategy)

solve_poker_strategy()