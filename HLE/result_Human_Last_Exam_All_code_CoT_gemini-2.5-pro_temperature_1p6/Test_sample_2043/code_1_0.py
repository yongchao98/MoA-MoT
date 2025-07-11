def solve_ems_scenario():
    """
    This script determines the best hospital destination based on a patient's
    critical condition by evaluating time and facility capabilities.
    """

    # Define the available hospital choices as a list of dictionaries
    choices = [
        {'id': 'A', 'name': 'Level 4 trauma center', 'time': 6, 'trauma_level': 4, 'can_handle_traumatic_arrest': False},
        {'id': 'B', 'name': 'Level 3 trauma center', 'time': 7, 'trauma_level': 3, 'can_handle_traumatic_arrest': True},
        {'id': 'C', 'name': 'Level 2 trauma center', 'time': 8, 'trauma_level': 2, 'can_handle_traumatic_arrest': True},
        {'id': 'D', 'name': 'Level 2 trauma center with a toxicologist', 'time': 15, 'trauma_level': 2, 'can_handle_traumatic_arrest': True},
        {'id': 'E', 'name': 'Level 1 trauma center with toxicologist', 'time': 15, 'trauma_level': 1, 'can_handle_traumatic_arrest': True},
    ]

    # Patient is in traumatic cardiac arrest. Time and surgical capability are paramount.
    # Set a maximum acceptable transport time for a patient in cardiac arrest.
    max_transport_time = 10  # minutes

    print("Step 1: Filtering out destinations that are too far away.")
    print(f"The patient is in cardiac arrest, so we eliminate any choice over {max_transport_time} minutes.")
    
    # Filter for time
    time_appropriate_choices = []
    for choice in choices:
        if choice['time'] <= max_transport_time:
            time_appropriate_choices.append(choice)
            print(f"- {choice['id']}: {choice['name']} ({choice['time']} minutes) is within the acceptable time frame.")
        else:
            print(f"- {choice['id']}: {choice['name']} ({choice['time']} minutes) is eliminated due to excessive transport time.")
    
    print("\nStep 2: Filtering the remaining options by capability.")
    print("A patient in traumatic cardiac arrest requires a facility with immediate surgical capabilities.")

    # Filter for capability
    fully_appropriate_choices = []
    for choice in time_appropriate_choices:
        if choice['can_handle_traumatic_arrest']:
            fully_appropriate_choices.append(choice)
            print(f"- {choice['id']}: {choice['name']} (Level {choice['trauma_level']}) is an appropriate facility.")
        else:
            print(f"- {choice['id']}: {choice['name']} (Level {choice['trauma_level']}) is eliminated because it cannot provide definitive surgical care.")

    print("\nStep 3: Selecting the best option from the suitable choices.")
    print("From the remaining appropriate facilities, choose the one with the shortest transport time.")

    # Find the choice with the minimum time among the fully appropriate ones
    best_choice = min(fully_appropriate_choices, key=lambda x: x['time'])

    print("\n--- Conclusion ---")
    print(f"The patient is in traumatic cardiac arrest. The primary need is immediate surgical intervention to address the cause of arrest.")
    print(f"Options D (15 minutes) and E (15 minutes) are too far.")
    print(f"Option A (6 minutes) is a Level 4 trauma center and is not equipped for the required surgery.")
    print(f"Between Option B (Level 3 at 7 minutes) and Option C (Level 2 at 8 minutes), both are capable.")
    print(f"Because every minute counts in cardiac arrest, the closest appropriate facility is the best choice.")
    print(f"Final Decision: {best_choice['id']} is the best destination as it is the closest hospital (at {best_choice['time']} minutes) with the necessary surgical capabilities (Level {best_choice['trauma_level']}).")
    
solve_ems_scenario()
<<<B>>>