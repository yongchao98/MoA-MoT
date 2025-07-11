import math

def solve_tournament():
    """
    Calculates the minimum number of days to find a winner in a tournament
    with 128 warriors, considering travel and fighting constraints.
    """
    initial_warriors = 128
    warriors = initial_warriors
    days = 0
    round_num = 1

    print("Tournament Plan:")
    print(f"Starting with {initial_warriors} warriors, each in their own city.\n")

    while warriors > 1:
        fights_this_round = warriors // 2
        
        print(f"--- Round {round_num} ({warriors} warriors -> {fights_this_round} winners) ---")
        
        # Day for Travel
        days += 1
        print(f"Day {days}: {fights_this_round} warriors travel to their opponents' cities.")
        
        # Day for Fighting
        days += 1
        print(f"Day {days}: {fights_this_round} fights take place in parallel. {fights_this_round} winners remain.")
        
        warriors = fights_this_round
        round_num += 1
        print("-" * 35 + "\n")

    num_rounds = math.log2(initial_warriors)
    days_per_round = 2
    total_days = num_rounds * days_per_round
    
    print("The winner has been determined!\n")
    print("Summary of Calculation:")
    print(f"Total number of rounds = log2({initial_warriors}) = {int(num_rounds)}")
    print(f"Days required per round = 1 (for travel) + 1 (for fighting) = {days_per_round}")
    print("Final equation for total minimum days:")
    print(f"{int(num_rounds)} * {days_per_round} = {int(total_days)}")

solve_tournament()
<<<14>>>