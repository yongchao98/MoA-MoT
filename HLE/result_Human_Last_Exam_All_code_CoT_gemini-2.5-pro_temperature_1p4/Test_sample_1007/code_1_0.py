def find_unchanged_ballet_step():
    """
    Analyzes classical ballet steps to find which one ends in the same
    leg position as it started.

    The analysis is stored in a dictionary where the key is the step name
    and the value is a boolean indicating if the feet change position
    (True = changes, False = remains the same).
    """

    # Data representing whether a step changes the starting foot position.
    # True: The foot that was in front is now in the back.
    # False: The foot that was in front is still in front at the end.
    steps_analysis = {
        "A. Entrechat six": False,
        "B. Échappé battu changé": True,
        "C. Assemblé": True,
        "D. Glissade derrière": True,
        "E. Gargouillade": True
    }

    print("Analyzing which ballet step starts and ends in the same position...")
    print("-" * 60)
    
    correct_step = None
    
    # Iterate through the steps to find the one that does not change position.
    for step, changes_position in steps_analysis.items():
        start_position = "Right Foot Front"
        # Determine the end position based on our analysis.
        if changes_position:
            end_position = "Left Foot Front"
        else:
            end_position = "Right Foot Front"

        print(f"Step: {step}")
        print(f"  Start Position: {start_position}")
        print(f"  End Position:   {end_position}")

        if start_position == end_position:
            correct_step = step
            print("  Result: Position is the SAME.\n")
        else:
            print("  Result: Position has CHANGED.\n")
            
    print("-" * 60)
    if correct_step:
        print(f"The correct answer is: {correct_step}")
    else:
        print("No step was found that meets the criteria.")

# Execute the function
find_unchanged_ballet_step()