import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to find a winner in a tournament
    with 128 warriors, starting in 128 different cities.
    """
    
    warriors = 128
    
    # 1. Determine the number of rounds in a single-elimination tournament.
    #    This is log base 2 of the number of warriors.
    num_rounds = int(math.log2(warriors))
    
    # 2. Determine the minimum days required for each round.
    #    - Day 1: Travel. To pair up N warriors, N/2 of them must travel to meet
    #      the other N/2. This takes 1 day.
    #    - Day 2: Fight. The N/2 matches can happen in parallel in N/2 cities.
    #      This takes 1 day.
    #    A warrior cannot travel and fight on the same day, and the winners of
    #    one round must be known before the next round's travel can begin.
    #    So, each round takes 2 days.
    days_per_round = 2
    
    # 3. Calculate the total minimum days.
    total_days = num_rounds * days_per_round
    
    # 4. Print the step-by-step explanation and the final equation.
    print(f"A tournament starts with {warriors} warriors in {warriors} different cities.")
    print("To determine a single winner, a single-elimination format is required.")
    print(f"\nNumber of elimination rounds needed = log2({warriors}) = {num_rounds} rounds.")
    print("\nEach round requires a minimum of 2 days:")
    print(f"  - 1 day for travel (so opponents can meet in the same city).")
    print(f"  - 1 day for fighting.")
    print("\nSince the rounds are sequential, the total time is the product of the number of rounds and the days per round.")
    
    print("\nFinal Equation:")
    print(f"Total Minimum Days = Number of Rounds * Days per Round")
    print(f"Total Minimum Days = {num_rounds} * {days_per_round} = {total_days}")

solve_tournament_days()