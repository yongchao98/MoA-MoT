import math

def calculate_tournament_duration():
    """
    Calculates the minimum number of days to determine a winner in a tournament
    with specific travel and fighting rules.
    """
    num_warriors = 128

    # A single-elimination tournament with N warriors requires log2(N) rounds for a winner to emerge.
    num_rounds = int(math.log2(num_warriors))

    # For any match to occur, at least one warrior must travel to the opponent's city.
    # Travel takes 1 day. A warrior can either travel or fight on any given day.
    # Let's analyze the time taken for each round from a participant's perspective.
    # Round 1:
    # Day 1: One warrior of a pair travels.
    # Day 2: The fight happens.
    # So, the first round is completed on Day 2 for the winners.
    # Each subsequent round also needs one day for travel and one day for fighting.
    days_per_round = 2
    
    # The total minimum number of days is the number of rounds multiplied by the days required per round.
    total_days = num_rounds * days_per_round

    print("--- Tournament Schedule Analysis ---")
    print(f"Number of Warriors: {num_warriors}")
    print(f"Number of Rounds (log2({num_warriors})): {num_rounds}")
    print(f"Minimum Days per Round (1 for travel + 1 for fighting): {days_per_round}")
    print("\n--- Final Calculation ---")
    print("The minimum number of days is the product of the number of rounds and the days required for each round.")
    print(f"Total Minimum Days = Number of Rounds * Days per Round")
    print(f"Total Minimum Days = {num_rounds} * {days_per_round}")
    print(f"Total Minimum Days = {total_days}")
    print("\nTherefore, the winner will be determined in a minimum of 14 days.")

calculate_tournament_duration()