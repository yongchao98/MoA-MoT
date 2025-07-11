import math

def solve_tournament_days():
    """
    Calculates the minimum number of days for a tournament with specific rules.
    """
    num_warriors = 128
    
    print(f"The tournament starts with {num_warriors} warriors, each in a different city.")
    print("The goal is to find the winner in the minimum number of days.\n")
    print("The optimal strategy involves parallel elimination rounds.")

    # The number of rounds in a single-elimination tournament is log2(N).
    if num_warriors <= 0 or (num_warriors & (num_warriors - 1) != 0 and num_warriors != 1):
        print("This calculation assumes the number of warriors is a power of 2 for a clean bracket.")
        return
        
    num_rounds = int(math.log2(num_warriors))

    print(f"A single-elimination tournament for {num_warriors} warriors will have {num_rounds} rounds.\n")
    print("Each round requires one day for travel and one day for fighting:")

    day_counter = 0
    current_warriors = num_warriors

    for i in range(num_rounds):
        round_num = i + 1
        
        # Travel Day
        day_counter += 1
        warriors_traveling = current_warriors // 2
        print(f"Day {day_counter}: (Round {round_num} Travel) {warriors_traveling} warriors travel to meet their opponents.")

        # Fight Day
        day_counter += 1
        fights_this_round = current_warriors // 2
        current_warriors //= 2
        print(f"Day {day_counter}: (Round {round_num} Fight) {fights_this_round} fights occur, leaving {current_warriors} warriors.\n")

    print("--- Summary ---")
    travel_days = num_rounds
    fight_days = num_rounds
    total_days = travel_days + fight_days

    print(f"Total number of rounds: {num_rounds}")
    print(f"Total days for travel: {travel_days} (1 day per round)")
    print(f"Total days for fighting: {fight_days} (1 day per round)")
    
    # Final equation as requested
    print("\nFinal calculation:")
    print(f"Total Minimum Days = {travel_days} + {fight_days} = {total_days}")

solve_tournament_days()