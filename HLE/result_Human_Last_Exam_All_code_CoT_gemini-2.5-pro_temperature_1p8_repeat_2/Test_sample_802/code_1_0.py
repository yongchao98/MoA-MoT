import math

def calculate_tournament_duration():
    """
    Calculates the minimum number of days to determine a winner in a tournament
    with 128 warriors, each starting in their own city.
    """
    
    num_warriors = 128
    days_per_round = 2  # 1 day for travel, 1 day for fighting

    # The tournament is single-elimination. The number of rounds required is
    # the logarithm base 2 of the number of warriors.
    # 128 -> 64 -> 32 -> 16 -> 8 -> 4 -> 2 -> 1 winner.
    num_rounds = math.log2(num_warriors)

    # The total minimum days is the number of rounds multiplied by the days per round.
    total_days = num_rounds * days_per_round
    
    print("To find the winner from 128 warriors, we need to complete a series of elimination rounds.")
    print("Each round consists of two phases that take one day each:")
    print("1. Travel Day: Half the warriors travel to their opponent's city.")
    print("2. Fight Day: Fights occur in parallel in different cities.")
    print("-" * 30)

    # We calculate the total number of rounds needed.
    print(f"Number of Warriors: {num_warriors}")
    print(f"Number of Rounds = log2({num_warriors}) = {int(num_rounds)}")
    
    # We calculate the total days.
    print(f"Days per Round = {days_per_round} (1 Travel + 1 Fight)")
    print("-" * 30)
    print("The final equation for the total minimum number of days is:")
    print(f"Total Days = Number of Rounds * Days per Round")
    print(f"Total Days = {int(num_rounds)} * {days_per_round} = {int(total_days)}")

calculate_tournament_duration()