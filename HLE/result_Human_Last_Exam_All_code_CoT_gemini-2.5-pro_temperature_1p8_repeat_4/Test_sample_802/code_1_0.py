import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to find a winner in a tournament.
    """
    
    # 1. Set the initial number of warriors.
    num_warriors = 128

    # 2. Determine the number of rounds in a single-elimination tournament.
    #    For N warriors, this is log base 2 of N.
    num_rounds = math.log2(num_warriors)
    if not num_rounds.is_integer():
        # This handles cases where the number of warriors is not a power of 2,
        # though for this problem it is.
        num_rounds = math.ceil(num_rounds)
    num_rounds = int(num_rounds)
    
    # 3. Determine the number of days required per round.
    #    Day 1: Travel. To have k/2 battles, k/2 warriors travel to their opponents' cities.
    #    Day 2: Fight. The k/2 battles occur in parallel in k/2 different cities.
    days_per_round = 2
    
    # 4. Calculate the total minimum days.
    total_days = num_rounds * days_per_round
    
    print("Step-by-step calculation:")
    print(f"Initial number of warriors: {num_warriors}")
    print(f"A single-elimination tournament requires log2({num_warriors}) rounds.")
    print(f"Number of rounds = {num_rounds}")
    print(f"Each round requires {days_per_round} days (1 for travel, 1 for fighting).")
    print(f"Final equation: Total days = {num_rounds} * {days_per_round}")
    print(f"The minimum number of days to determine the winner is: {total_days}")

solve_tournament_days()