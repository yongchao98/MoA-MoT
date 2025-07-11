def find_max_diameter_cutoff():
    """
    Analyzes clinical data to find the maximum Aortomesenteric diameter cutoff
    that meets specific sensitivity and specificity criteria for identifying
    esophagogastric varices.
    """
    # Step 1: Define the data from the clinical study.
    # Each dictionary represents a cutoff value and its diagnostic performance.
    # The 'cutoff_mm' is interpreted as "diameter <= value".
    data = [
        {'cutoff_mm': 9.8, 'sensitivity': 72.5, 'specificity': 72.1},
        {'cutoff_mm': 8.5, 'sensitivity': 60.9, 'specificity': 83.7},
        {'cutoff_mm': 7.5, 'sensitivity': 47.8, 'specificity': 93.0}
    ]

    # Step 2: Define the required performance criteria.
    required_sensitivity = 60
    required_specificity = 80

    print("Searching for the maximum Aortomesenteric diameter cutoff with:")
    print(f"Sensitivity > {required_sensitivity}%")
    print(f"Specificity > {required_specificity}%")
    print("-" * 30)

    # Step 3: Evaluate each cutoff and find the ones that are valid.
    valid_cutoffs = []
    print("Evaluating available data points:")
    for item in data:
        is_sensitive_enough = item['sensitivity'] > required_sensitivity
        is_specific_enough = item['specificity'] > required_specificity
        
        print(f"\nTesting cutoff <= {item['cutoff_mm']} mm:")
        print(f"  - Sensitivity check: {item['sensitivity']}% > {required_sensitivity}% -> {is_sensitive_enough}")
        print(f"  - Specificity check: {item['specificity']}% > {required_specificity}% -> {is_specific_enough}")

        if is_sensitive_enough and is_specific_enough:
            valid_cutoffs.append(item['cutoff_mm'])
            print("  - Result: This cutoff is VALID.")
        else:
            print("  - Result: This cutoff is INVALID.")

    # Step 4: Find the maximum valid cutoff from the list.
    if not valid_cutoffs:
        print("\nNo cutoff value satisfies both criteria.")
        final_answer = None
    else:
        max_valid_cutoff = max(valid_cutoffs)
        print("-" * 30)
        print(f"The valid cutoff value(s) are: {valid_cutoffs}")
        print(f"The maximum valid cutoff value is {max_valid_cutoff}.")
        print("This means the test 'diameter <= 8.5 mm' meets the criteria.")
        print(f"Therefore, the maximum Aortomesenteric diameter that will be identified by this test is {max_valid_cutoff} mm.")
        final_answer = max_valid_cutoff

    # Final explicit output of the number in the equation context
    if final_answer is not None:
        print(f"\nThe final equation can be stated as: Maximum qualifying diameter = {final_answer}")

if __name__ == '__main__':
    find_max_diameter_cutoff()