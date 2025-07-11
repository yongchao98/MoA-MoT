import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to find a winner in the tournament.
    """
    num_warriors = 128

    # --- Plan Explanation ---
    print("To find the minimum number of days for the tournament, we need the most efficient schedule.")
    print("The strategy is to maximize the number of fights happening in parallel.")
    print("\nStep 1: Determine the number of elimination rounds.")
    print("The tournament starts with 128 warriors and ends with 1. Each round halves the number of warriors.")
    # Number of rounds = log base 2 of the number of warriors
    # 2^x = 128 => x = 7
    num_rounds = math.log2(num_warriors)
    print(f"The number of rounds required is log2({num_warriors}) = {int(num_rounds)} rounds.")

    print("\nStep 2: Determine the time taken for each round.")
    days_per_round = 2  # 1 day for travel + 1 day for fighting
    print("Each round consists of two phases:")
    print("  - Travel (1 day): Half the warriors travel to meet their opponents.")
    print("  - Fighting (1 day): The paired-up warriors fight in parallel.")
    print(f"So, each round takes {days_per_round} days.")

    print("\nStep 3: Calculate the total minimum days.")
    # Total days = (Number of rounds) * (Days per round)
    total_days = num_rounds * days_per_round
    
    # Final equation output as requested
    print("The final calculation is:")
    print(f"Total Minimum Days = Number of Rounds * Days per Round")
    print(f"Total Minimum Days = {int(num_rounds)} * {days_per_round}")
    print(f"Total Minimum Days = {int(total_days)}")

solve_tournament_days()