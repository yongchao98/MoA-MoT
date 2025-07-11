def find_max_diameter_threshold():
    """
    Finds the maximum Aortomesenteric diameter threshold that is greater than 
    60% sensitive and greater than 80% specific for identifying EVP enhancement.
    
    This function uses a representative dataset for demonstration.
    """
    # Step 1: Establish the Data
    # Data represents: (Diameter Threshold (mm), Sensitivity (%), Specificity (%))
    # This data is hypothetical but represents a realistic diagnostic trade-off.
    diagnostic_data = [
        (4.0, 95, 30),
        (5.0, 90, 45),
        (6.0, 85, 60),
        (7.0, 75, 75),
        (8.0, 68, 82),
        (9.0, 62, 88),
        (10.0, 55, 95),
        (11.0, 45, 98)
    ]

    # Step 2: Define Conditions
    min_sensitivity = 60
    min_specificity = 80
    
    print(f"Searching for a diameter threshold with sensitivity > {min_sensitivity}% and specificity > {min_specificity}%\n")

    valid_thresholds = []
    
    # Step 3: Iterate and Filter
    for diameter, sensitivity, specificity in diagnostic_data:
        print(f"Checking diameter threshold > {diameter} mm...")
        
        is_sensitive = sensitivity > min_sensitivity
        is_specific = specificity > min_specificity
        
        print(f"  - Sensitivity check: {sensitivity} > {min_sensitivity} ({is_sensitive})")
        print(f"  - Specificity check: {specificity} > {min_specificity} ({is_specific})")

        if is_sensitive and is_specific:
            valid_thresholds.append(diameter)
            print(f"  -> Result: Condition met. {diameter} mm is a valid threshold.\n")
        else:
            print(f"  -> Result: Condition not met.\n")
            
    # Step 4: Find the Maximum
    if not valid_thresholds:
        final_answer = None
        print("No diameter threshold met the specified criteria.")
    else:
        max_valid_diameter = max(valid_thresholds)
        final_answer = max_valid_diameter
        
        # Find the corresponding sensitivity and specificity for the final answer
        final_sensitivity = 0
        final_specificity = 0
        for diameter, sensitivity, specificity in diagnostic_data:
            if diameter == final_answer:
                final_sensitivity = sensitivity
                final_specificity = specificity
                break
        
        # Step 5: Display the Result
        print("="*50)
        print("Final analysis:")
        print(f"The valid diameter thresholds are: {valid_thresholds}")
        print(f"The maximum valid threshold is {final_answer} mm.")
        print(f"At this threshold, the final equation check is:")
        print(f"Sensitivity: {final_sensitivity} > {min_sensitivity}")
        print(f"Specificity: {final_specificity} > {min_specificity}")

if __name__ == '__main__':
    find_max_diameter_threshold()