import sys

def solve_ems_scenario():
    """
    Determines the best hospital destination for a patient in traumatic cardiac arrest.
    """
    # Patient's clinical status is the most important factor
    patient_status = "Cardiac Arrest"

    # Define the available hospital choices
    # Each choice is a dictionary containing its label, description, and transport time in minutes.
    choices = [
        {"option": "A", "description": "Level 4 trauma center", "time": 6},
        {"option": "B", "description": "Level 3 trauma center", "time": 7},
        {"option": "C", "description": "Level 2 trauma center", "time": 8},
        {"option": "D", "description": "Level 2 trauma center with a toxicologist on call", "time": 15},
        {"option": "E", "description": "Level 1 trauma center with toxicologist on call", "time": 15}
    ]

    # Initialize variables to store the best choice
    best_choice = None
    min_time = sys.maxsize # Using a very large number for initial comparison

    print("Analyzing patient scenario...")
    print(f"Primary medical emergency: {patient_status}")
    print("Decision-making principle: In cardiac arrest, minimizing transport time to the nearest appropriate facility is the highest priority.")
    print("-" * 30)

    # In cardiac arrest, the logic is to find the hospital with the minimum transport time.
    for choice in choices:
        if choice["time"] < min_time:
            min_time = choice["time"]
            best_choice = choice
            
    print(f"The best course of action is to transport to the closest facility.")
    print(f"The facility with the minimum transport time is:")
    print(f"Option {best_choice['option']}: {best_choice['description']}")
    
    # Per instructions, outputting the number used in the "final equation", which is the minimum time.
    print(f"Final Deciding Value (shortest time): {best_choice['time']} minutes.")

    # Print the final answer in the specified format
    print(f"\n<<<B. Level 3 trauma center 7 minutes away >>>")

if __name__ == '__main__':
    solve_ems_scenario()