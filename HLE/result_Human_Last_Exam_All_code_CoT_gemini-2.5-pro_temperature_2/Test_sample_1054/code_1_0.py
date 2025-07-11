def find_max_diameter_threshold():
    """
    Finds the maximum Aortomesenteric (AM) diameter threshold that meets specific
    sensitivity and specificity criteria for identifying EVP enhancement.
    """
    # Step 1: Define the clinical data.
    # This data represents different AM diameter thresholds and their corresponding
    # sensitivity and specificity percentages.
    # Format: [AM Diameter (mm), Sensitivity (%), Specificity (%)]
    clinical_data = [
        [4.0, 85, 70],
        [4.5, 80, 81],
        [5.0, 75, 86],
        [5.5, 68, 90],
        [6.0, 61, 94],
        [6.5, 54, 97],
        [7.0, 45, 99]
    ]

    # Step 2: Define the required criteria.
    min_sensitivity = 60
    min_specificity = 80

    # Step 3: Identify all diameter thresholds that satisfy the criteria.
    valid_thresholds = []
    for diameter, sensitivity, specificity in clinical_data:
        if sensitivity > min_sensitivity and specificity > min_specificity:
            valid_thresholds.append([diameter, sensitivity, specificity])

    # Step 4: Find the threshold with the maximum diameter from the valid list.
    if not valid_thresholds:
        print("No diameter threshold was found that meets the specified criteria.")
        return

    # Find the entry with the maximum diameter among the valid ones
    max_valid_threshold_entry = max(valid_thresholds, key=lambda x: x[0])
    
    max_diameter = max_valid_threshold_entry[0]
    final_sensitivity = max_valid_threshold_entry[1]
    final_specificity = max_valid_threshold_entry[2]
    
    # Step 5: Print the results and the logic.
    print(f"Finding the maximum Aortomesenteric diameter that is > {min_sensitivity}% sensitive and > {min_specificity}% specific.")
    print("-" * 60)
    print("Valid thresholds found (Sensitivity > 60% and Specificity > 80%):")
    for d, s, sp in valid_thresholds:
        print(f"  - Diameter: {d} mm, Sensitivity: {s}%, Specificity: {sp}%")
    print("-" * 60)
    
    print("The maximum diameter from this valid list is chosen.")
    print("\nFinal Answer Derivation:")
    print(f"The Aortomesenteric diameter of {max_diameter} mm results in:")
    print(f"  - Sensitivity: {final_sensitivity}% (which is > {min_sensitivity}%)")
    print(f"  - Specificity: {final_specificity}% (which is > {min_specificity}%)")
    
    print(f"\nThis is the highest threshold that meets both conditions. Therefore, the answer is {max_diameter}.")

if __name__ == '__main__':
    find_max_diameter_threshold()