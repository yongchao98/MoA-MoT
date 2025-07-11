import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to complete a tournament
    with the given constraints.
    """
    num_warriors = 128

    # To have one winner, we need to eliminate (num_warriors - 1) warriors.
    # In a single-elimination tournament, the number of rounds is log base 2
    # of the number of participants.
    num_rounds = int(math.log2(num_warriors))

    # Each round requires two phases:
    # 1. Travel: One warrior from each pair must travel to their opponent's city.
    # This takes 1 day.
    # 2. Fighting: The matches can happen in parallel in different cities.
    # This takes another day, as a warrior cannot travel and fight on the same day.
    days_per_round_travel = 1
    days_per_round_fight = 1
    total_days_per_round = days_per_round_travel + days_per_round_fight

    # The total minimum number of days is the number of rounds multiplied by
    # the number of days required for each round.
    min_total_days = num_rounds * total_days_per_round

    print("To determine the winner from 128 warriors, we need a single-elimination tournament.")
    print(f"Since 128 = 2^7, the tournament will have {num_rounds} rounds.")
    print("\nEach round requires a minimum of 2 days:")
    print(f"- {days_per_round_travel} day for warriors to travel to their opponent's city.")
    print(f"- {days_per_round_fight} day for the actual fight.")
    
    print("\nThe final equation for the total minimum number of days is:")
    print("Number of Rounds * (Days for Travel + Days for Fighting) = Total Days")
    
    # Printing each number in the final equation as requested.
    print(f"{num_rounds} * ({days_per_round_travel} + {days_per_round_fight}) = {min_total_days}")

solve_tournament_days()