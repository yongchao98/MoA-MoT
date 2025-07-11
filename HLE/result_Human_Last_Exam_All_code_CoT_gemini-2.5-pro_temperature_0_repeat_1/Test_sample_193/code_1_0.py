def solve_poker_strategy():
    """
    Calculates the optimal poker strategy for the given scenario.
    """
    # 1. Define game variables
    pot = 10
    hero_value_combos = 1  # AA
    hero_bluff_combos_available = 10 # QQ-33

    # 2. Determine the optimal number of bluffs to use
    # To maximize EV, we must use the maximum number of available bluffs.
    # The strategy involves betting with all value hands and all chosen bluff hands.
    num_bluffs_to_bet = hero_bluff_combos_available

    # 3. Calculate the optimal bet size to make the villain indifferent
    # The formula is B = Pot * (Number of Bluffs / Number of Value Hands)
    bet_size = pot * (num_bluffs_to_bet / hero_value_combos)
    bet_size = int(bet_size) # The result is a whole number

    # 4. Determine the strategy for AA
    # The strategy that maximizes our EV involves always betting our value hands.
    # A deviation (checking AA) would result in a lower EV.
    # Therefore, the frequency of betting with AA is 100%.
    action = "BET"
    frequency = 100

    # Rounding to the nearest even number as per instructions
    if bet_size % 2 != 0:
        bet_size = round(bet_size / 2) * 2
    if frequency % 2 != 0:
        frequency = round(frequency / 2) * 2

    # 5. Print the results and the equation used
    print("This is a game theory problem where we must construct a betting range that is difficult for our opponent to play against.")
    print("Our goal is to choose a bet size and a number of bluffs that maximizes our profit.")
    print("The optimal strategy is to bet with all our value hands (AA) and all our available bluff hands (QQ-33).")
    print("\nThe bet size is calculated to make the villain indifferent between calling and folding.")
    print(f"The equation is: Bet Size = Pot * (Number of Bluffs / Number of Value Combos)")
    
    # Output each number in the final equation as requested
    print(f"Calculation: {bet_size} = {pot} * ({num_bluffs_to_bet} / {hero_value_combos})")
    
    print("\nSince we always want to take this highest EV action with our value hand AA, our strategy for AA is pure.")
    print("\nFinal Answer:")
    
    # The final answer format
    final_answer_string = f"{action} {bet_size} {frequency}%"
    print(final_answer_string)
    
    # Required for the platform
    print(f"\n<<<{final_answer_string}>>>")

solve_poker_strategy()