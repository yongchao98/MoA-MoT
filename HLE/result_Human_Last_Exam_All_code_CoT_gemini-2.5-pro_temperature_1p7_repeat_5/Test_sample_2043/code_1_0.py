import sys

def solve_ems_dilemma():
    """
    Analyzes hospital options for a critical patient and recommends the best destination.
    """
    # Patient status: Traumatic cardiac arrest (highest priority), Tylenol overdose (secondary).
    # This requires weighting Trauma Level highest.

    # Weights for scoring criteria
    trauma_level_weight = 10
    toxicology_support_weight = 2
    transport_time_penalty = -1 # Negative because more time is bad

    # Map trauma centers to a numerical score for calculation
    trauma_level_map = {
        1: 5,  # Highest capability
        2: 4,
        3: 3,
        4: 1,  # Lowest capability (stabilize and transfer)
    }

    options = {
        'A': {'name': 'Level 4 trauma center', 'level': 4, 'time': 6, 'has_tox': False},
        'B': {'name': 'Level 3 trauma center', 'level': 3, 'time': 7, 'has_tox': False},
        'C': {'name': 'Level 2 trauma center', 'level': 2, 'time': 8, 'has_tox': False},
        'D': {'name': 'Level 2 trauma center with a toxicologist', 'level': 2, 'time': 15, 'has_tox': True},
        'E': {'name': 'Level 1 trauma center with toxicologist', 'level': 1, 'time': 15, 'has_tox': True}
    }

    best_option = None
    max_score = -sys.maxsize # Start with a very low score

    print("Evaluating EMS destination options based on weighted scores:")
    print("Formula: (Trauma Level Score * {t_weight}) + (Transport Time * {time_pen}) + (Toxicology Support * {tox_weight})\n".format(
        t_weight=trauma_level_weight,
        time_pen=transport_time_penalty,
        tox_weight=toxicology_support_weight
    ))

    for key, attrs in options.items():
        trauma_score = trauma_level_map[attrs['level']]
        tox_score = 1 if attrs['has_tox'] else 0
        time_val = attrs['time']

        total_score = (trauma_score * trauma_level_weight) + \
                      (time_val * transport_time_penalty) + \
                      (tox_score * toxicology_support_weight)
        
        print("Option {key}: {name}".format(key=key, name=attrs['name']))
        # The prompt requires showing the final equation with numbers
        print("Score = ({t_score} * {t_weight}) + ({time} * {time_pen}) + ({tox} * {tox_weight}) = {total}".format(
            t_score=trauma_score, t_weight=trauma_level_weight,
            time=time_val, time_pen=transport_time_penalty,
            tox=tox_score, tox_weight=toxicology_support_weight,
            total=total_score
        ))
        print("-" * 20)

        if total_score > max_score:
            max_score = total_score
            best_option = key

    print("\nConclusion: The patient's most critical condition is traumatic cardiac arrest.")
    print("This requires the highest level of care, making Trauma Level the most important factor.")
    print("Option {best} has the highest score of {score}, making it the best destination despite the longer transport time.".format(
        best=best_option, score=max_score
    ))


solve_ems_dilemma()
<<<E>>>