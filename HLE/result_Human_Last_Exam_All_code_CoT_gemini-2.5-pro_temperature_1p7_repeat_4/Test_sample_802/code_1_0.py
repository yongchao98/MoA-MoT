import math

def solve_tournament_riddle():
    """
    Calculates the minimum number of days to determine a winner in a tournament
    with 128 warriors, based on the given rules.
    """
    initial_warriors = 128
    days_per_round = 2 # 1 day for travel, 1 day for fighting

    # 1. Determine the number of rounds.
    # The tournament is single-elimination. The number of rounds required to get
    # from N participants to 1 winner is log2(N).
    num_rounds = math.log2(initial_warriors)
    if not num_rounds.is_integer():
        # This case handles scenarios where the number of warriors is not a power of 2,
        # but for 128, it is an integer.
        num_rounds = math.ceil(num_rounds)
    num_rounds = int(num_rounds)
    
    # 2. Calculate the total days.
    # Each round takes 2 days: one for travel and one for fighting.
    total_days = num_rounds * days_per_round
    
    print("Step 1: Calculate the number of tournament rounds.")
    print(f"With {initial_warriors} warriors, the number of rounds is log2({initial_warriors}).")
    print(f"Number of rounds = {num_rounds}")
    print("-" * 30)

    print("Step 2: Calculate the days required per round.")
    print("Each round needs warriors to be in the same city to fight.")
    print("Day 1 (Travel): Half the warriors travel to meet the other half.")
    print("Day 2 (Fight): The paired warriors fight.")
    print(f"Days per round = {days_per_round}")
    print("-" * 30)
    
    print("Step 3: Calculate the total minimum days.")
    print("Total Days = (Number of Rounds) * (Days per Round)")
    # The final equation as requested
    print(f"The final equation is: {num_rounds} * {days_per_round} = {total_days}")
    print("-" * 30)
    
    print(f"The minimum number of days to determine the winner is {total_days}.")

solve_tournament_riddle()
<<<14>>>