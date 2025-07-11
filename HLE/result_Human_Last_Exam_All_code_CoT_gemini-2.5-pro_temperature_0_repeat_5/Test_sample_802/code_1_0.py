import math

def calculate_tournament_duration():
    """
    Calculates and explains the minimum number of days for the tournament.
    """
    initial_warriors = 128
    
    print(f"Calculating the minimum days for a tournament with {initial_warriors} warriors.")
    print("Each warrior starts in a different city, and travel between any two cities takes 1 day.")
    print("A warrior can either travel or fight on a given day.\n")

    # Calculate the number of rounds
    # Number of rounds = log2(number of warriors)
    num_rounds = int(math.log2(initial_warriors))
    
    print(f"The tournament is a single-elimination format, which will require {num_rounds} rounds.\n")
    print("--- Tournament Schedule Breakdown ---")

    current_warriors = initial_warriors
    total_days = 0
    
    # Loop through each round to explain the process
    for i in range(1, num_rounds + 1):
        matches_in_round = current_warriors // 2
        
        # Day for Travel
        total_days += 1
        print(f"Round {i}:")
        print(f"  Day {total_days}: {matches_in_round} warriors travel to their opponents' cities.")
        
        # Day for Fighting
        total_days += 1
        print(f"  Day {total_days}: {matches_in_round} fights occur, eliminating {matches_in_round} warriors.")
        
        current_warriors //= 2
        print(f"  {current_warriors} winners remain.\n")

    print("--- Final Calculation ---")
    days_per_round = 2
    final_days = num_rounds * days_per_round
    print("Each of the 7 rounds requires 1 day for travel and 1 day for fighting.")
    print(f"Total minimum days = {num_rounds} rounds * {days_per_round} days/round = {final_days} days.")

calculate_tournament_duration()