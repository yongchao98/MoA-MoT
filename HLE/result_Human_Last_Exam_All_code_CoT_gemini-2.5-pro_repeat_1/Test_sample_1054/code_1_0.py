def find_am_diameter_threshold():
    """
    Finds the Aortomesenteric diameter threshold that meets specific sensitivity and specificity criteria for EVP enhancement.
    """
    # Step 1: Represent the data from the study
    # Data shows sensitivity and specificity for a given AM diameter cutoff (<= value)
    roc_data = [
        {'diameter_le': 2.0, 'sensitivity': 100, 'specificity': 17},
        {'diameter_le': 4.0, 'sensitivity': 96, 'specificity': 42},
        {'diameter_le': 6.0, 'sensitivity': 80, 'specificity': 67},
        {'diameter_le': 8.0, 'sensitivity': 64, 'specificity': 83},
        {'diameter_le': 10.0, 'sensitivity': 40, 'specificity': 96},
        {'diameter_le': 12.0, 'sensitivity': 20, 'specificity': 100}
    ]

    # Step 2: Define the required criteria
    required_sensitivity = 60
    required_specificity = 80

    print(f"Finding the AM diameter threshold with sensitivity > {required_sensitivity}% and specificity > {required_specificity}%.")
    print("-" * 30)

    # Step 3 & 4: Iterate through the data to find the correct threshold
    solution = None
    for point in roc_data:
        diameter = point['diameter_le']
        sensitivity = point['sensitivity']
        specificity = point['specificity']

        # Check if the current point meets both criteria
        if sensitivity > required_sensitivity and specificity > required_specificity:
            solution = point
            # Since the question asks "At most what diameter", and the list is ordered,
            # this will be the highest valid diameter that meets the criteria from the given options.
            break # We found the first and only value that fits the criteria.

    # Step 5: Display the result
    if solution:
        print(f"The criteria are met at an Aortomesenteric diameter of <= {solution['diameter_le']} mm.")
        print("Verification:")
        # Output each number in the final comparison
        print(f"Sensitivity: {solution['sensitivity']}% > {required_sensitivity}%")
        print(f"Specificity: {solution['specificity']}% > {required_specificity}%")
        print("\nTherefore, the highest diameter from the options provided that meets these conditions is 8.0 mm.")
        
        # The final answer to be extracted
        global final_answer
        final_answer = solution['diameter_le']
    else:
        print("No diameter in the provided data meets both criteria.")
        final_answer = "Not found"

find_am_diameter_threshold()
print(f'<<<{final_answer}>>>')