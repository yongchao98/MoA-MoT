import math

def calculate_minimum_tournament_days():
    """
    Calculates and explains the minimum number of days to determine a winner
    in a tournament of 128 warriors with specific travel and battle constraints.
    """
    num_warriors = 128

    # 1. Determine the number of rounds for a single-elimination tournament.
    # For N warriors, the number of rounds is log2(N).
    num_rounds = int(math.log2(num_warriors))

    # 2. Determine the number of days required per round.
    # Each round requires pairing up the winners from the previous round.
    # Since all warriors start in different cities, and winners of a round
    # will also be in different cities, pairing requires travel.
    travel_days_per_round = 1  # Half the warriors travel to meet the other half.
    fight_days_per_round = 1   # Battles happen the next day, in parallel.
    
    # Each round takes (travel + fight) days because a warrior cannot
    # do both on the same day.
    days_per_round = travel_days_per_round + fight_days_per_round

    # 3. Calculate the total minimum days.
    total_days = num_rounds * days_per_round

    # 4. Print the explanation and the final result.
    print("Step-by-step calculation for the minimum tournament duration:")
    print(f"  - Initial number of warriors: {num_warriors}")
    
    print("\nStep 1: Calculate the number of tournament rounds.")
    print(f"  - The tournament is single-elimination, so the number of rounds is log2({num_warriors}).")
    print(f"  - Number of rounds = {num_rounds}")

    print("\nStep 2: Calculate the days required per round.")
    print(f"  - For each round, winners must be paired. One warrior must travel to the other's city.")
    print(f"  - Travel time = {travel_days_per_round} day.")
    print(f"  - The next day, the battles can occur.")
    print(f"  - Battle time = {fight_days_per_round} day.")
    print(f"  - Total days per round = {travel_days_per_round} + {fight_days_per_round} = {days_per_round} days.")
    
    print("\nStep 3: Calculate the total minimum days for the tournament.")
    print(f"  - Total days = (Number of rounds) * (Days per round)")
    print(f"  - Final Equation: {num_rounds} * {days_per_round} = {total_days}")

calculate_minimum_tournament_days()