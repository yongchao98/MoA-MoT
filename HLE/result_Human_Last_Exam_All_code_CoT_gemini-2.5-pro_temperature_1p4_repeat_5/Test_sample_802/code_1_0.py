import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to find a winner in a tournament.
    """
    num_warriors = 128
    
    # To get a single winner from N warriors, we need N-1 eliminations.
    # The most efficient way is a single-elimination tournament structure,
    # which proceeds in rounds.

    # 1. Calculate the number of rounds.
    # Each round halves the number of warriors.
    # The number of rounds is log base 2 of the number of warriors.
    num_rounds = math.log2(num_warriors)
    
    # 2. Calculate the days per round.
    # For each round, half of the remaining warriors must travel to meet their opponents.
    travel_days_per_round = 1
    # On the next day, the fights can happen in parallel.
    fight_days_per_round = 1
    days_per_round = travel_days_per_round + fight_days_per_round
    
    # 3. Calculate the total minimum days.
    # Total days = (Number of rounds) * (Days per round)
    total_days = num_rounds * days_per_round
    
    print("Step 1: Determine the number of tournament rounds.")
    print(f"With {int(num_warriors)} warriors, we need to find the number of rounds (R) in a single-elimination bracket.")
    print(f"R = log2({int(num_warriors)}) = {int(num_rounds)} rounds\n")
    
    print("Step 2: Determine the minimum days per round.")
    print("For each round, half of the warriors must travel to their opponents' cities.")
    print(f"Travel takes {travel_days_per_round} day.")
    print("The day after travel, the fights can occur.")
    print(f"Fighting takes {fight_days_per_round} day.")
    print(f"Total days per round = {travel_days_per_round} (travel) + {fight_days_per_round} (fight) = {days_per_round} days\n")
    
    print("Step 3: Calculate the total minimum days for the tournament.")
    print(f"Total days = Number of rounds * Days per round")
    print(f"Total days = {int(num_rounds)} * {days_per_round} = {int(total_days)}")

solve_tournament_days()
<<<14>>>