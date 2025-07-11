import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to find a winner in a tournament.
    """
    # 1. Define the initial number of warriors.
    num_warriors = 128

    # 2. Calculate the number of rounds in a single-elimination tournament.
    # This is log base 2 of the number of warriors.
    num_rounds = math.log2(num_warriors)
    if not num_rounds.is_integer():
        # This case is for when the number of warriors is not a power of 2,
        # but the problem specifies 128.
        num_rounds = math.ceil(num_rounds)
    num_rounds = int(num_rounds)

    # 3. Define the days required for travel and fighting in each round.
    # To start a round, half the warriors must travel to meet their opponents.
    days_for_travel = 1
    # The next day, the fights can happen in parallel in different cities.
    days_for_fighting = 1
    
    # Each round consists of a travel day and a fight day.
    days_per_round = days_for_travel + days_for_fighting

    # 4. Calculate the total minimum days for the tournament.
    total_minimum_days = num_rounds * days_per_round

    # 5. Print the step-by-step calculation.
    print(f"The tournament starts with {num_warriors} warriors.")
    print(f"The number of elimination rounds required is log2({num_warriors}), which is {num_rounds} rounds.")
    print(f"For each round, one day is needed for travel ({days_for_travel}) and one day is needed for fighting ({days_for_fighting}).")
    print(f"This means each round takes {days_for_travel} + {days_for_fighting} = {days_per_round} days.")
    print("\nFinal Calculation:")
    print(f"Total Minimum Days = {num_rounds} (rounds) * {days_per_round} (days per round)")
    print(f"Total Minimum Days = {total_minimum_days}")

solve_tournament_days()