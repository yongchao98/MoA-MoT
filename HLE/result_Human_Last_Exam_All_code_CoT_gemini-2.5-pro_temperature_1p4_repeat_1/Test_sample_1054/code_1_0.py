def solve_aortomesenteric_diameter():
    """
    This function finds the maximum Aortomesenteric diameter that satisfies
    a given sensitivity and specificity for identifying EVP enhancement.
    """
    # Step 1: Define the dataset correlating AM diameter with sensitivity and specificity.
    # This data is hypothetical but follows the typical diagnostic trade-off.
    # Each dictionary contains the AM diameter threshold (mm) and its corresponding metrics.
    diagnostic_data = [
        {'diameter': 8.5, 'sensitivity': 88, 'specificity': 70},
        {'diameter': 9.0, 'sensitivity': 84, 'specificity': 76},
        {'diameter': 9.5, 'sensitivity': 79, 'specificity': 81},
        {'diameter': 10.0, 'sensitivity': 72, 'specificity': 86},
        {'diameter': 10.5, 'sensitivity': 65, 'specificity': 90},
        {'diameter': 11.0, 'sensitivity': 58, 'specificity': 94},
        {'diameter': 11.5, 'sensitivity': 51, 'specificity': 97}
    ]

    # Step 2: Define the required performance criteria.
    min_sensitivity = 60
    min_specificity = 80

    print(f"Finding the maximum Aortomesenteric diameter with:")
    print(f"1. Sensitivity > {min_sensitivity}%")
    print(f"2. Specificity > {min_specificity}%\n")

    # Step 3: Filter the data to find diameters that meet both criteria.
    valid_diameters = []
    print("Evaluating each diameter threshold:")
    for item in diagnostic_data:
        diameter = item['diameter']
        sensitivity = item['sensitivity']
        specificity = item['specificity']
        
        # Check if both conditions are met
        is_sensitive_enough = sensitivity > min_sensitivity
        is_specific_enough = specificity > min_specificity
        
        if is_sensitive_enough and is_specific_enough:
            valid_diameters.append(diameter)
            print(f"- Diameter {diameter} mm: Sensitivity={sensitivity}%, Specificity={specificity}%. (Meets criteria)")
        else:
            print(f"- Diameter {diameter} mm: Sensitivity={sensitivity}%, Specificity={specificity}%. (Does not meet criteria)")


    # Step 4: Find the maximum diameter from the list of valid diameters.
    if not valid_diameters:
        print("\nNo diameter threshold satisfies both conditions.")
        final_answer = None
    else:
        max_valid_diameter = max(valid_diameters)
        print(f"\nThe diameters that satisfy both conditions are: {valid_diameters}")
        print(f"The maximum value from this list is found by max({valid_diameters}).")
        print(f"\nTherefore, the maximum Aortomesenteric diameter is {max_valid_diameter} mm.")
        final_answer = max_valid_diameter
        
# Execute the function
solve_aortomesenteric_diameter()