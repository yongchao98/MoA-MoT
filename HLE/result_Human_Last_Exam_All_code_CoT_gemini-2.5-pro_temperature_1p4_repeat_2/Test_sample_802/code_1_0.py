import math

def calculate_tournament_days():
    """
    Calculates the minimum number of days to determine a winner in a tournament
    with 128 warriors, considering travel and fight constraints.
    """
    initial_warriors = 128
    warriors = initial_warriors
    days = 0
    round_count = 1

    print("My plan is to simulate the tournament round by round to find the total number of days.")
    print(f"The tournament starts with {warriors} warriors, each in a separate city.\n")

    # A round consists of 1 day for travel and 1 day for fighting.
    while warriors > 1:
        fights = warriors // 2
        print(f"--- Round {round_count} ({warriors} -> {warriors // 2} warriors) ---")

        # Day for travel: half the warriors travel to meet their opponents.
        days += 1
        print(f"Day {days}: {fights} warriors travel to their opponents' cities.")
        
        # Day for fighting: fights occur in parallel in different cities.
        days += 1
        print(f"Day {days}: {fights} fights happen, eliminating {fights} warriors.")
        
        warriors //= 2
        print(f"End of Round {round_count}: {warriors} warriors remain. Total days elapsed: {days}.\n")
        round_count += 1
        
    print("--- Final Calculation ---")
    num_rounds = math.log2(initial_warriors)
    days_per_round = 2
    print(f"The tournament needs {int(num_rounds)} rounds to go from {initial_warriors} warriors to 1.")
    print(f"Each round takes {days_per_round} days (1 for travel + 1 for fighting).")
    print("\nThe final equation is:")
    print(f"{int(num_rounds)} * {days_per_round} = {days}")
    print(f"\nThe minimum number of days to determine the winner is {days}.")

calculate_tournament_days()