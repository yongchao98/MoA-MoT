import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to complete a tournament
    with 128 warriors under the specified rules.
    """
    
    warriors = 128

    # 1. Calculate the number of rounds in a single-elimination tournament.
    # For N participants, where N is a power of 2, this is log2(N).
    num_rounds = int(math.log2(warriors))

    # 2. Calculate the number of days required per round.
    # Each round requires bringing pairs of warriors together. Since they start
    # in different cities, this involves travel.
    # Day X: One warrior from each pair travels to meet the other. (1 day)
    # Day X+1: The matches for the round take place in parallel. (1 day)
    days_per_round = 2 # 1 day for travel + 1 day for fighting

    # 3. Calculate the total minimum days for the tournament.
    total_days = num_rounds * days_per_round

    print("The tournament starts with 128 warriors, so it will have a total of 7 rounds (128 -> 64 -> 32 -> 16 -> 8 -> 4 -> 2 -> 1).")
    print("\nEach round requires two days to complete:")
    print(" - 1 Day for Travel: Half the remaining warriors travel to meet their opponents.")
    print(" - 1 Day for Fighting: All matches of the round occur in parallel in different cities.")
    print(f"\nThus, each of the {num_rounds} rounds takes {days_per_round} days.")
    
    print("\nThe total minimum number of days is calculated as:")
    print(f"{num_rounds} rounds * {days_per_round} days/round = {total_days} days")

solve_tournament_days()