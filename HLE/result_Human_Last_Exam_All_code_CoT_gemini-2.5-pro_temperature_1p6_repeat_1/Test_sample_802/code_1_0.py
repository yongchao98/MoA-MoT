import math

def solve_tournament_days():
    """
    Calculates the minimum number of days to determine a winner in a tournament.
    """
    num_warriors = 128

    # The tournament is a single-elimination format. To find a single winner from
    # N warriors, we need N-1 fights. The most efficient way to schedule these
    # is in rounds.

    # Step 1: Calculate the number of rounds.
    # In each round, the number of warriors is halved.
    # The number of rounds is the logarithm base 2 of the number of warriors.
    # For N warriors, it takes log2(N) rounds.
    num_rounds = int(math.log2(num_warriors))

    # Step 2: Calculate the days per round.
    # At the start of each round, the remaining warriors are in different cities.
    # To conduct a round of fights, we must pair them up.
    # - Day 1 (Travel): For each pair, one warrior travels to the other's city.
    #   Since this can happen for all pairs at the same time, this phase takes 1 day.
    # - Day 2 (Fight): The fights happen in the cities where the pairs met.
    #   Since each city has an arena, all fights can happen in parallel. This takes 1 day.
    # A warrior cannot travel and fight on the same day, so a round takes 2 days.
    days_per_round = 2

    # Step 3: Calculate the total minimum days.
    # Total Days = (Number of Rounds) * (Days per Round)
    total_days = num_rounds * days_per_round

    print(f"Number of warriors: {num_warriors}")
    print(f"Number of rounds needed: log2({num_warriors}) = {num_rounds}")
    print(f"Days required per round: 1 (for travel) + 1 (for fighting) = {days_per_round}")
    print("---")
    print("Final Calculation:")
    print(f"Total Minimum Days = {num_rounds} rounds * {days_per_round} days/round")
    print(f"Result: {total_days}")

solve_tournament_days()
<<<14>>>