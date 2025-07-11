import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to determine a winner in a tournament.
    """
    num_warriors = 128
    days_per_round = 2 # 1 day for travel, 1 day for fighting

    # The number of rounds in a single-elimination tournament is log2(number of participants).
    num_rounds = int(math.log2(num_warriors))

    # Total days = number of rounds * days per round.
    total_days = num_rounds * days_per_round

    print(f"The tournament starts with {num_warriors} warriors.")
    print("To find a single winner, the tournament is structured in rounds where half the warriors are eliminated each round.")
    print(f"The number of rounds required is log2({num_warriors}) = {num_rounds} rounds.")
    print("\nEach round requires a minimum of 2 days:")
    print("1. One day for half the warriors to travel to their opponents' cities.")
    print("2. One day for the fights to occur in parallel.")
    print("\nTherefore, the total minimum number of days to find the winner is:")
    print(f"{num_rounds} rounds * {days_per_round} days/round = {total_days} days")

solve_tournament_days()