import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to determine a winner in a tournament
    with 128 warriors, considering travel and fight constraints.
    """
    initial_warriors = 128

    # The tournament is single-elimination. The number of rounds required is
    # the number of times you can halve the participants until one remains.
    # This is calculated with log base 2.
    num_rounds = int(math.log2(initial_warriors))

    # Each round requires two phases because a warrior can't travel and fight on the same day.
    # 1. Travel Phase: To pair up for N/2 fights, N/2 warriors must travel. This takes 1 day.
    days_for_travel_per_round = 1
    # 2. Fight Phase: The N/2 fights can happen in parallel, each in a different city. This takes 1 day.
    days_for_fights_per_round = 1
    
    days_per_round = days_for_travel_per_round + days_for_fights_per_round

    # Total time is the product of the number of rounds and the time each round takes.
    total_days = num_rounds * days_per_round

    print(f"Initial number of warriors: {initial_warriors}")
    print("To determine a single winner, the tournament must proceed in elimination rounds.")
    print("-" * 30)
    
    print(f"Step 1: Calculate the number of rounds.")
    print(f"The number of rounds is log2({initial_warriors}) = {num_rounds} rounds.")
    print("\nStep 2: Calculate the days per round.")
    print(f"Each round requires half the warriors to travel (takes {days_for_travel_per_round} day) and then fight (takes {days_for_fights_per_round} day).")
    print(f"So, each round takes {days_for_travel_per_round} + {days_for_fights_per_round} = {days_per_round} days.")
    print("-" * 30)

    print("Step 3: Calculate the total minimum days.")
    print("Total Days = Number of Rounds * Days per Round")
    
    # Final equation with the numbers included
    print(f"Total Days = {num_rounds} * {days_per_round} = {total_days}")
    print("-" * 30)

    print(f"The winner will be determined in a minimum of {total_days} days.")

solve_tournament_days()