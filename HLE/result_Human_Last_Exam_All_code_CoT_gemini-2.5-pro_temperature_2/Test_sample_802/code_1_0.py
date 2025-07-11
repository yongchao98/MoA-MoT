import math

def calculate_tournament_duration():
    """
    Calculates the minimum number of days to determine a winner in the tournament.
    """
    
    # Total number of warriors at the start of the tournament.
    num_warriors = 128
    
    # Each round of a single-elimination tournament halves the number of participants.
    # The number of rounds is the logarithm base 2 of the number of warriors.
    # 2^7 = 128, so there are 7 rounds.
    num_rounds = int(math.log2(num_warriors))
    
    # Each round requires two phases: travel and fighting.
    # 1. Travel Phase: Half of the warriors travel to their opponents' cities. This takes 1 day.
    # 2. Fight Phase: All matches of the round occur in parallel in different cities. This takes 1 day.
    # Thus, each round takes 2 days.
    days_per_round = 2
    
    # The total minimum number of days is the number of rounds multiplied by the days per round.
    total_days = num_rounds * days_per_round
    
    # Print the final equation with the numbers plugged in.
    print(f"{num_rounds} * {days_per_round} = {total_days}")

calculate_tournament_duration()