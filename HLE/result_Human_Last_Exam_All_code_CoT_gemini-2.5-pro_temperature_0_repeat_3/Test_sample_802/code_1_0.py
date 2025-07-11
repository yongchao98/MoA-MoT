import math

def solve_tournament_days():
    """
    Calculates and explains the minimum number of days to determine a winner
    in a tournament of 128 warriors.
    """
    initial_warriors = 128
    warriors = initial_warriors
    days = 0
    round_num = 1

    print(f"A tournament begins with {initial_warriors} warriors, each in a separate city.")
    print("The optimal strategy is to pair warriors up and have them fight in parallel across different cities.")
    print("-" * 40)

    # Loop through the rounds until a single winner remains
    while warriors > 1:
        print(f"Round {round_num}: Starting with {warriors} warriors.")

        # Travel Day: Half the warriors travel to meet their opponents.
        days += 1
        travelers = warriors // 2
        print(f"  Day {days}: {travelers} warriors travel for one day to meet their opponents.")

        # Fight Day: Fights happen in parallel.
        days += 1
        fights = warriors // 2
        warriors_before_fight = warriors
        warriors = warriors // 2
        print(f"  Day {days}: {fights} fights occur in parallel, eliminating {fights} warriors.")
        print(f"    {warriors_before_fight} / 2 = {warriors} winners remain.")
        print("-" * 40)

        round_num += 1

    print(f"After {round_num - 1} rounds, a single winner is determined.")
    print(f"The total minimum number of days required is {days}.")

    # Final summary equation
    print("\n--- Final Calculation ---")
    num_rounds = math.log2(initial_warriors)
    days_per_round = 2  # 1 day for travel + 1 day for fighting
    total_days = num_rounds * days_per_round
    print(f"Number of Rounds = log2({initial_warriors}) = {int(num_rounds)}")
    print(f"Days per Round = 1 (Travel) + 1 (Fight) = {days_per_round}")
    print(f"Total Days = Number of Rounds * Days per Round")
    print(f"Total Days = {int(num_rounds)} * {days_per_round} = {int(total_days)}")

solve_tournament_days()