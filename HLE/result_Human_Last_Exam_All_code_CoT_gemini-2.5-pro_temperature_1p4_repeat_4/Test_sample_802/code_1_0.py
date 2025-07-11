import math

def solve_tournament():
    """
    Calculates the minimum number of days to find a winner in a tournament
    with 128 warriors, considering travel and battle constraints.
    """
    initial_warriors = 128
    current_warriors = initial_warriors
    total_days = 0
    round_number = 0
    
    print(f"Tournament starting with {initial_warriors} warriors.")
    print("Each round requires one day for travel and one day for battles.\n")

    # The list to hold the days per round for the final equation
    days_per_round_list = []

    while current_warriors > 1:
        round_number += 1
        days_this_round = 2 # 1 for travel, 1 for battle
        
        print(f"--- Round {round_number} ---")
        print(f"Warriors at start: {current_warriors}")
        
        # Day for Travel
        travelers = current_warriors // 2
        print(f"Day {total_days + 1}: {travelers} warriors travel to other cities.")
        
        # Day for Battles
        battles = current_warriors // 2
        print(f"Day {total_days + 2}: {battles} battles occur in parallel.")

        current_warriors //= 2
        total_days += days_this_round
        days_per_round_list.append(str(days_this_round))
        print(f"Warriors remaining: {current_warriors}\n")
        
    print("--- Tournament Complete ---")
    print("The total number of days is the sum of the days for each round.")

    # Building and printing the final equation
    equation_str = " + ".join(days_per_round_list)
    print(f"Final Calculation: {equation_str} = {total_days}")

solve_tournament()
<<<14>>>