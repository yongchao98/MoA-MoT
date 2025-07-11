import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to determine a winner in a tournament
    with 128 warriors, considering travel and fighting logistics.
    """
    num_warriors = 128

    print(f"A tournament starts with {num_warriors} warriors, each in their own city.")
    print("To find a single winner, a single-elimination tournament is required.")
    
    # Step 1: Calculate the number of rounds in the tournament.
    # In each round, the number of warriors is halved. The number of rounds is log2(number of warriors).
    num_rounds = math.log2(num_warriors)
    print(f"\nStep 1: Calculate the number of rounds.")
    print(f"   Number of rounds = log2({num_warriors}) = {int(num_rounds)} rounds.")

    # Step 2: Determine the minimum number of days per round.
    # To conduct a round of fights, pairs of warriors must be in the same city.
    # This requires one warrior from each pair to travel.
    travel_days_per_round = 1
    # Once the pairs are together, the fights can occur. Since all fights in a round
    # can happen in parallel in different cities, this takes one day.
    fight_days_per_round = 1
    days_per_round = travel_days_per_round + fight_days_per_round
    
    print(f"\nStep 2: Calculate the minimum days per round.")
    print(f"   Each round requires travel for half the warriors ({travel_days_per_round} day) and then a day for the fights ({fight_days_per_round} day).")
    print(f"   Minimum days per round = {travel_days_per_round} + {fight_days_per_round} = {days_per_round} days.")

    # Step 3: Calculate the total minimum days for the tournament.
    total_days = num_rounds * days_per_round
    print(f"\nStep 3: Calculate the total minimum days for the tournament.")
    print("   Total minimum days = (Number of rounds) * (Days per round).")

    print("\nFinal Calculation:")
    print(f"   {int(num_rounds)} * ({travel_days_per_round} + {fight_days_per_round}) = {int(total_days)}")
    print(f"\nThe winner will be determined in a minimum of {int(total_days)} days.")

solve_tournament_days()