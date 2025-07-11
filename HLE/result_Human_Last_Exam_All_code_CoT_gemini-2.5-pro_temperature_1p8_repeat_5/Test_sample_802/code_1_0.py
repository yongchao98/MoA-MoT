import math

def solve_tournament():
    """
    Calculates the minimum number of days to find a winner in a tournament
    with 128 warriors, considering travel and fighting constraints.
    """
    num_warriors = 128

    # --- Plan ---
    # 1. A single-elimination tournament reduces the number of participants by half each round.
    #    The number of rounds is log base 2 of the number of warriors.
    # 2. Each round requires two phases:
    #    a) Travel: To pair up for a fight, warriors must be in the same city. Since all
    #       warriors start in different cities, half of them must travel to meet the
    #       other half. This takes 1 day.
    #    b) Fight: Once the pairs are in the same cities, the fights can occur. Since each
    #       pair is in a different city with its own arena, all fights for a round can
    #       happen on the same day. This takes 1 day.
    # 3. The total time is the number of rounds multiplied by the days per round (1 travel + 1 fight).

    # Step 1: Calculate the number of rounds.
    num_rounds = math.log2(num_warriors)

    # Step 2: Calculate the days per round.
    days_for_travel = 1
    days_for_fighting = 1
    days_per_round = days_for_travel + days_for_fighting

    # Step 3: Calculate the total minimum days.
    total_days = num_rounds * days_per_round

    print("Step-by-step calculation for the tournament:")
    print("-" * 40)
    print(f"Initial number of warriors: {num_warriors}")
    print(f"The tournament is single-elimination, requiring {num_warriors - 1} fights in total.")
    print("-" * 40)
    
    print(f"1. Number of rounds required: log2({num_warriors}) = {int(num_rounds)} rounds.")
    
    print("\n2. Days required per round:")
    print(f"   - To pair up, half the warriors must travel: {days_for_travel} day.")
    print(f"   - After travel, fights happen in parallel: {days_for_fighting} day.")
    print(f"   - Total days per round = {days_for_travel} + {days_for_fighting} = {days_per_round} days.")

    print("-" * 40)
    print("3. Total minimum days for the tournament:")
    print("   Total Days = (Number of Rounds) * (Days per Round)")
    print(f"   Final Equation: {int(num_rounds)} * {days_per_round} = {int(total_days)}")
    print("-" * 40)
    print(f"\nThe minimum number of days to determine the winner is {int(total_days)}.")

solve_tournament()