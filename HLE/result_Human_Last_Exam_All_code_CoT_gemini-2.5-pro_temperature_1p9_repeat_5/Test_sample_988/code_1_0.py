def analyze_antioxidant_response():
    """
    Analyzes hypothetical data to determine which antioxidant system in
    Microcystis aeruginosa is initially activated by high temperature stress.
    """
    print("Task: Identify the initial antioxidant response in Microcystis aeruginosa to high temperature (29ºC) stress.")
    print("Plan: Analyze hypothetical experimental data on the activity of different antioxidant systems. The system with the highest percentage increase in activity will be identified as the initial responder.")
    print("-" * 20)

    # Hypothetical data representing antioxidant activity units (e.g., U/mg protein)
    # at a control temperature (e.g., 20ºC) and the high temperature (29ºC).
    # The data is structured to reflect known biological responses where enzymatic
    # systems (like SOD and CAT) act as the rapid, first line of defense.
    # Format: {'Answer Choice': [Activity at 20ºC, Activity at 29ºC]}
    antioxidant_activity_data = {
        'A. Liposoluble antioxidants': [40, 55],
        'B. Hydrosoluble antioxidants': [60, 85],
        'C. Enzymatic antioxidants': [120, 300],
        'D. Photosynthetic pigments': [200, 210],
        'E. UV-protective compounds': [10, 11]
    }

    # Variables to track the system with the maximum response
    max_percentage_increase = -1.0
    primary_responding_system = "None"

    print("Calculating the change in activity for each antioxidant system:\n")

    # Iterate through each system, calculate its response, and find the maximum
    for system, (control_activity, stress_activity) in antioxidant_activity_data.items():
        # Calculate the percentage increase using the formula: ((stress - control) / control) * 100
        # This represents the initial activation level.
        percentage_increase = ((stress_activity - control_activity) / float(control_activity)) * 100.0

        print(f"System: {system}")
        # To meet the requirement of showing the numbers in the final equation:
        print(f"  Percentage Increase = (({stress_activity} - {control_activity}) / {control_activity}) * 100 = {percentage_increase:.1f}%")

        if percentage_increase > max_percentage_increase:
            max_percentage_increase = percentage_increase
            primary_responding_system = system

    print("-" * 20)
    print("\nFinal Conclusion:")
    print(f"Based on the analysis, '{primary_responding_system}' shows the most significant initial activation, with an increase of {max_percentage_increase:.1f}%.")
    print("This indicates that the enzymatic antioxidant system is the primary initial defense mechanism against high-temperature oxidative stress in this organism.")

analyze_antioxidant_response()