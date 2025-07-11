import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to complete a tournament
    with 128 warriors, explaining each step of the calculation.
    """
    num_warriors = 128
    days_per_round = 2 # 1 day for travel, 1 day for fighting

    # 1. Calculate the number of rounds
    # In a single-elimination tournament, the number of participants is halved each round.
    # The number of rounds is log2 of the number of warriors.
    num_rounds = math.log2(num_warriors)

    # 2. Calculate total days
    # Total days = (Number of Rounds) * (Days per Round)
    total_days = num_rounds * days_per_round

    print("Step 1: Determine the number of tournament rounds.")
    print(f"A tournament with {num_warriors} warriors requires log2({num_warriors}) rounds.")
    print(f"Number of rounds = {int(num_rounds)}")
    print("-" * 30)

    print("Step 2: Determine the minimum days per round.")
    print("For each round, the winners from the previous round are in different cities.")
    print("To fight, one warrior from each pair must travel to their opponent's city.")
    print("Day 1: Travel (half the warriors travel).")
    print("Day 2: Fight (all pairs fight in parallel in different cities).")
    print(f"Minimum days per round = {days_per_round}")
    print("-" * 30)

    print("Step 3: Calculate the total minimum days for the tournament.")
    print("Total Days = Number of Rounds * Days per Round")
    # The final equation with each number explicitly shown
    print(f"Total Days = {int(num_rounds)} * {days_per_round} = {int(total_days)}")

solve_tournament_days()