import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to determine a winner in a tournament
    with 128 warriors, considering travel and fighting constraints.
    """
    
    # The initial number of warriors
    num_warriors = 128
    
    # In a single-elimination tournament, the number of participants is halved each round.
    # The total number of rounds is the logarithm base 2 of the number of warriors.
    num_rounds = math.log2(num_warriors)
    
    # For each round, we analyze the time required.
    # Let's take the first round (128 warriors -> 64 winners) as an example.
    # We need to conduct 64 fights. The warriors start in 128 different cities.
    
    # Step 1: Travel. To form 64 pairs in 64 cities, 64 warriors must travel to meet
    # their opponents. This can happen in parallel and takes 1 day.
    travel_days_per_round = 1
    
    # Step 2: Fight. The 64 fights can happen in parallel on the next day, each in a
    # separate city. This takes 1 day. A warrior cannot travel and fight on the same day.
    fight_days_per_round = 1
    
    # The total number of days per round is the sum of travel and fight days.
    days_per_round = travel_days_per_round + fight_days_per_round
    
    # The total minimum time is the number of rounds multiplied by the days per round.
    total_days = num_rounds * days_per_round
    
    print("### Tournament Logistics Breakdown ###")
    print(f"Initial number of warriors: {num_warriors}")
    
    print("\nStep 1: Determine the number of tournament rounds.")
    print(f"The number of rounds is log2({num_warriors}) = {int(num_rounds)} rounds.")
    
    print("\nStep 2: Determine the number of days per round.")
    print(f"Each round requires a travel day and a fight day.")
    print(f"Days per round = {travel_days_per_round} (travel) + {fight_days_per_round} (fight) = {days_per_round} days.")
    
    print("\nStep 3: Calculate the total minimum days for the tournament.")
    print("Total Days = Number of Rounds * Days per Round")
    
    # The final equation with all numbers, as requested.
    print(f"The calculation is: {int(num_rounds)} * {days_per_round} = {int(total_days)}")
    
solve_tournament_days()