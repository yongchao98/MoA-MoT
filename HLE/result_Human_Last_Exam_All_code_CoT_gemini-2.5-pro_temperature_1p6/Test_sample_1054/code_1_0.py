import sys

def find_max_diameter_cutoff():
    """
    Finds the maximum Aortomesenteric diameter cutoff based on sensitivity and specificity criteria.
    The underlying medical condition is assumed to be Superior Mesenteric Artery (SMA) syndrome.
    """
    # Step 1: Data from a clinical study on Aortomesenteric Diameter for diagnosing SMA syndrome.
    # A measurement <= the cutoff 'diameter' is considered a positive test.
    # Data is adapted from relevant literature (e.g., Ozkurt et al., 2007).
    # Format: {cutoff_diameter_mm: (sensitivity_%, specificity_%)}
    clinical_data = {
        6.0: (42.8, 100.0),
        7.0: (57.1, 100.0),
        8.0: (100.0, 100.0),
        9.0: (100.0, 85.7),
        10.0: (100.0, 57.1),
        11.0: (100.0, 50.0),
        12.0: (100.0, 42.8)
    }

    # Step 2: Define the required performance criteria.
    required_sensitivity = 60.0
    required_specificity = 80.0

    print("Finding the maximum Aortomesenteric diameter cutoff with:")
    print(f"Sensitivity > {required_sensitivity}% AND Specificity > {required_specificity}%")
    print("-" * 50)

    valid_diameters = []

    # Step 3 & 4: Iterate through the data, check each cutoff, and print the logic.
    sorted_diameters = sorted(clinical_data.keys())
    for diameter in sorted_diameters:
        sensitivity, specificity = clinical_data[diameter]
        
        # Check if the cutoff meets the criteria
        is_sensitive_enough = sensitivity > required_sensitivity
        is_specific_enough = specificity > required_specificity

        print(f"Checking cutoff <= {diameter} mm:")
        # "Equation" part: explicitly showing the numbers and the boolean check
        print(f"  Condition 1: Sensitivity ({sensitivity}%) > {required_sensitivity}%  -->  Result: {is_sensitive_enough}")
        print(f"  Condition 2: Specificity ({specificity}%) > {required_specificity}%  -->  Result: {is_specific_enough}")
        
        if is_sensitive_enough and is_specific_enough:
            valid_diameters.append(diameter)
            print("  Conclusion: Cutoff is VALID.\n")
        else:
            print("  Conclusion: Cutoff is INVALID.\n")

    # Step 5: Determine the maximum valid diameter from the list.
    if not valid_diameters:
        print("No diameter cutoff meets the specified criteria.")
        max_diameter = None
    else:
        max_diameter = max(valid_diameters)
        print("-" * 50)
        print(f"The set of valid diameter cutoffs is: {valid_diameters} mm")
        print(f"The maximum value from this set is the answer.")
        print(f"\nAt most, the Aortomesenteric diameter cutoff can be {max_diameter} mm to satisfy the criteria.")

    if max_diameter is not None:
        sys.stdout.write(f'<<<{max_diameter}>>>')

if __name__ == '__main__':
    find_max_diameter_cutoff()