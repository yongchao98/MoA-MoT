def find_worst_opening_hand():
    """
    Analyzes a poker scenario to find the worst suited Jack to open from the button.
    """
    scenario = "100BB deep, Rake-Free Cash Game, Open-raise from Button to 2.5x"
    position = "Button (BTN)"
    hand_category = "Suited Jacks (Jxs)"

    print(f"Analyzing the poker hand range for the following scenario:")
    print(f"Scenario: {scenario}")
    print("-" * 40)

    print("Step 1: Evaluating the Button's opening strategy.")
    print(f"In poker, the {position} is the most powerful pre-flop position.")
    print("This allows a player to profitably open-raise with a very wide range of hands.")
    print("-" * 40)

    print("Step 2: Considering the 'Rake-Free' condition.")
    print("A rake-free game increases the profitability of marginal hands.")
    print("Therefore, the optimal opening range is even wider than in a standard game.")
    print("-" * 40)

    print("Step 3: Determining the GTO (Game Theory Optimal) range.")
    print(f"For this scenario, the GTO opening range from the {position} is very wide.")
    print(f"This range includes all {hand_category}.")
    all_suited_jacks = ["JAs", "JKs", "JQs", "JTs", "J9s", "J8s", "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"]
    print(f"The suited Jack hands in this opening range are: {', '.join(all_suited_jacks)}.")
    print("-" * 40)

    print("Step 4: Identifying the 'worst' hand in the opening range.")
    print("The 'worst' hand is the one with the lowest kicker, as it has the least equity.")
    worst_hand = all_suited_jacks[-1] # The last one in the sorted list
    print(f"Out of all the suited Jacks that should be opened, the worst one is {worst_hand}.")
    print("-" * 40)

    print(f"The final answer is: {worst_hand}")


if __name__ == "__main__":
    find_worst_opening_hand()
