import math

def solve_tournament():
    """
    Calculates the minimum number of days to complete a tournament
    with 128 warriors, explaining the logic step-by-step.
    """
    initial_warriors = 128
    warriors = initial_warriors
    total_days = 0
    round_num = 1
    
    # This list will store the '2' for each round to form the final equation
    days_per_round_list = []

    print(f"Starting a tournament with {initial_warriors} warriors, each in their own city.\n")

    while warriors > 1:
        fights_in_round = warriors // 2
        
        # Each round requires 1 day for travel and 1 day for fights
        days_for_this_round = 2
        
        print(f"--- Round {round_num} ---")
        print(f"Warriors at start of round: {warriors}")
        print(f"Matches to be played: {fights_in_round}")
        print(f"Day {total_days + 1}: {fights_in_round} warriors travel to meet their opponents.")
        print(f"Day {total_days + 2}: {fights_in_round} fights occur in parallel in {fights_in_round} cities.")
        
        # Update state for the next round
        warriors //= 2
        total_days += days_for_this_round
        days_per_round_list.append(str(days_for_this_round))
        
        print(f"Warriors remaining: {warriors}")
        print(f"Total days elapsed: {total_days}\n")
        round_num += 1
        
    num_rounds = round_num - 1
    
    print("--- Tournament Complete ---")
    print(f"To determine a single winner, {num_rounds} rounds of elimination were needed.")
    print("Each round took 2 days (1 for travel, 1 for fighting).")
    print("\nThe final calculation for the total minimum days is:")
    
    # Creates the equation string like "2 + 2 + ... + 2"
    equation = " + ".join(days_per_round_list)
    
    print(f"{equation} = {total_days}")

solve_tournament()