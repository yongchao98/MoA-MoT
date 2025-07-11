import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to complete a tournament
    with 128 warriors, explaining the logic step-by-step.
    """
    initial_warriors = 128
    warriors = initial_warriors
    days = 0
    round_num = 1

    print(f"A tournament with {initial_warriors} warriors begins.")
    print("The strategy is to conduct the tournament in knockout rounds to maximize parallel fights.")
    print("-" * 40)

    # Loop through each round of the tournament until a winner is found
    while warriors > 1:
        print(f"Round {round_num}:")
        print(f"  - Warriors at start of round: {int(warriors)}")
        
        fights_in_round = warriors / 2
        
        # Day for travel
        days += 1
        print(f"  - Day {days}: {int(fights_in_round)} warriors travel for one day to meet their opponents.")
        
        # Day for fighting
        days += 1
        print(f"  - Day {days}: {int(fights_in_round)} fights happen in parallel across {int(fights_in_round)} cities.")
        
        # Update the number of remaining warriors
        warriors /= 2
        print(f"  - Warriors remaining: {int(warriors)}")
        print("-" * 40)
        round_num += 1
        
    print("The final winner has been determined.")
    print(f"Total minimum number of days required: {days}")
    print("\nSummary Calculation:")
    
    # Calculate the number of rounds using log base 2
    num_rounds = math.log2(initial_warriors)
    days_per_round = 2
    
    print(f"Number of rounds = log2({initial_warriors}) = {int(num_rounds)}")
    print(f"Days per round = 1 (for travel) + 1 (for fighting) = {days_per_round}")
    
    # Display the final equation with all numbers
    print("\nFinal Equation:")
    print(f"{int(num_rounds)} rounds * {days_per_round} days/round = {int(num_rounds * days_per_round)} days")

solve_tournament_days()