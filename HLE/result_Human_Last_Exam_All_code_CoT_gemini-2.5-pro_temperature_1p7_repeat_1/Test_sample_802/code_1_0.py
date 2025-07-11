import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to determine a winner in a tournament
    with 128 warriors, considering travel and fighting constraints.
    """
    num_warriors = 128
    
    # A single-elimination tournament with N warriors has log2(N) rounds.
    if num_warriors <= 0 or (num_warriors & (num_warriors - 1) != 0):
        print("The number of warriors must be a power of 2 for this calculation.")
        return

    num_rounds = int(math.log2(num_warriors))
    
    # Each round requires pairing up warriors. Since they start in different cities,
    # this involves travel followed by a fight.
    # Day 1 of a round: Half the warriors travel to meet their opponents.
    # Day 2 of a round: The matches take place.
    # Therefore, each round takes 2 days.
    days_per_round = 2
    
    total_days = num_rounds * days_per_round
    
    print(f"Starting a tournament with {num_warriors} warriors.")
    print("The tournament is single-elimination, requiring a series of rounds.")
    print("-" * 50)
    
    print("Calculation Breakdown:")
    # The prompt requires printing each number in the final equation.
    print(f"1. Determine the number of rounds:")
    print(f"   Number of Rounds = log2({num_warriors}) = {num_rounds} rounds")
    
    print("\n2. Determine the days needed per round:")
    print(f"   Each round requires 1 day for travel and 1 day for fighting.")
    print(f"   Days per Round = 1 (travel) + 1 (fight) = {days_per_round} days")
    
    print("\n3. Calculate the total minimum days:")
    print(f"   Total Days = (Number of Rounds) * (Days per Round)")
    print(f"   Total Days = {num_rounds} * {days_per_round} = {total_days}")
    print("-" * 50)
    
    print(f"The minimum number of days to determine the winner is {total_days}.")

solve_tournament_days()