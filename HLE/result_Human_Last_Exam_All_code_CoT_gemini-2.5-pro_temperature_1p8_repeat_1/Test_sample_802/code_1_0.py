import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to determine a winner in a tournament.
    """
    # Initial number of warriors, each in their own city.
    initial_warriors = 128

    print("Step-by-step breakdown of the tournament:")
    print("=" * 45)

    current_warriors = initial_warriors
    day_counter = 0
    round_number = 1

    # The tournament proceeds in rounds until only one warrior is left.
    while current_warriors > 1:
        print(f"Round {round_number}: Starting with {current_warriors} warriors.")
        
        # In each round, half the warriors must travel to meet the other half.
        # This takes one full day.
        travelers = current_warriors // 2
        day_counter += 1
        print(f"  - Day {day_counter}: {travelers} warriors travel to meet their opponents.")
        
        # After traveling, the pairs of warriors are in the same cities.
        # The matches can now happen in parallel, each in a different city.
        # This takes another full day.
        matches = current_warriors // 2
        day_counter += 1
        print(f"  - Day {day_counter}: {matches} matches are held in {matches} different cities.")
        
        # After the matches, half the warriors are eliminated.
        current_warriors //= 2
        print(f"  => End of Round {round_number}: {current_warriors} warriors remain.")
        print("-" * 45)
        round_number += 1
    
    print("The winner has been determined.")
    print("\nFinal Calculation:")
    
    # The number of rounds is log2 of the number of warriors.
    num_rounds = math.log2(initial_warriors)
    # Each round consists of 1 travel day and 1 fight day.
    days_per_round = 2
    
    # Total days = Number of Rounds * Days per Round
    total_days = num_rounds * days_per_round
    
    print("The minimum number of days is the number of rounds multiplied by the days per round (1 for travel + 1 for fighting).")
    print(f"Total Days = log2({initial_warriors}) * {days_per_round} = {int(num_rounds)} * {days_per_round} = {int(total_days)}")

solve_tournament_days()