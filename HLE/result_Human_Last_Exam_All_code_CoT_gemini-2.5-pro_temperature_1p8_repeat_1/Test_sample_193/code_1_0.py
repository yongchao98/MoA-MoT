import math

# Define the core parameters of the poker scenario
pot_size = 10
stack_size = 1000
value_combos = 6 # Number of AA combos
bluff_candidates = 60 # Number of QQ-33 combos

def solve_poker_strategy():
    """
    Calculates and explains the optimal GTO strategy for AA in the given scenario.
    """
    # In this scenario, our hand (AA) is the nuts and beats the villain's hand (KK).
    # The goal is to devise a strategy for AA that maximizes our Expected Value (EV).

    # The optimal strategy involves betting as large as possible to maximize value.
    # The maximum possible bet is the effective stack size.
    optimal_bet_size = stack_size

    # This large bet must be balanced with bluffs to be unexploitable.
    # We calculate the number of bluff combos needed to make the villain indifferent
    # to calling our all-in bet.
    # The formula is: Num_Bluffs = Num_Value_Combos * Bet / (Pot + Bet)
    bluffs_needed = value_combos * optimal_bet_size / (pot_size + optimal_bet_size)

    # Since this betting strategy maximizes the EV of our nut hand (AA), it is the
    # optimal play. We should therefore do it 100% of the time.
    bet_frequency = 100

    # Per the instructions, round the final sizing and percentage to the nearest even number.
    rounded_bet_size = int(round(optimal_bet_size / 2.0) * 2)
    rounded_frequency = int(round(bet_frequency / 2.0) * 2)

    # Now, we print the explanation and the final answer.
    print("Poker Strategy Analysis:")
    print("=" * 25)
    print(f"Pot: ${pot_size}, Stacks: ${stack_size}")
    print("My Hand (AA) vs. Villain's Hand (KK) on a 2s2d2c board.")
    print("-" * 25)
    print("Goal: Maximize profit (EV) with our nut hand, AA.")
    print("Logic: The Expected Value of a value bet increases with the bet size.")
    print(f"Optimal Strategy: Bet the maximum possible amount, which is all-in for ${stack_size}.")
    print("-" * 25)
    print("Verifying the Strategy (Balancing Equation):")
    print("Num_Bluffs_Needed = Num_Value_Combos * Bet_Size / (Pot_Size + Bet_Size)")
    print(f"Num_Bluffs_Needed = {value_combos} * {optimal_bet_size} / ({pot_size} + {optimal_bet_size}) = {bluffs_needed:.2f}")
    print(f"This is a valid strategy as we have {bluff_candidates} combos available for bluffing.")
    print("-" * 25)

    final_strategy_string = f"BET {rounded_bet_size} {rounded_frequency}%"
    print("Final Optimal Strategy for AA:")
    print(final_strategy_string)

# Execute the function to solve the task
solve_poker_strategy()