import sys

def find_max_diameter_threshold():
    """
    This function analyzes a set of diagnostic accuracy data to find the maximum
    Aortomesenteric (AM) diameter threshold that meets specific sensitivity and
    specificity criteria for identifying EVP enhancement.
    
    The criteria are:
    - Sensitivity must be greater than 60%
    - Specificity must be greater than 80%
    """
    
    # Hypothetical data from a medical study.
    # The 'threshold' indicates that a measurement less than or equal to
    # this value is considered a positive test for the condition.
    roc_data = [
        {'threshold': 4.0, 'sensitivity': 0.25, 'specificity': 0.99},
        {'threshold': 5.0, 'sensitivity': 0.45, 'specificity': 0.95},
        {'threshold': 6.0, 'sensitivity': 0.58, 'specificity': 0.91},
        {'threshold': 7.0, 'sensitivity': 0.72, 'specificity': 0.88},
        {'threshold': 8.0, 'sensitivity': 0.85, 'specificity': 0.75},
        {'threshold': 9.0, 'sensitivity': 0.95, 'specificity': 0.60},
    ]

    # Required minimum values for sensitivity and specificity
    min_sensitivity = 0.60
    min_specificity = 0.80

    valid_thresholds = []

    print("Analyzing thresholds...\n")
    # Iterate through each data point to find valid thresholds
    for item in roc_data:
        threshold = item['threshold']
        sensitivity = item['sensitivity']
        specificity = item['specificity']

        # Check if the current threshold meets the criteria
        is_sensitivity_valid = sensitivity > min_sensitivity
        is_specificity_valid = specificity > min_specificity
        
        print(f"Checking threshold <= {threshold} mm:")
        # Print the check for sensitivity
        print(f"  Is sensitivity {sensitivity} > {min_sensitivity}? {is_sensitivity_valid}")
        # Print the check for specificity
        print(f"  Is specificity {specificity} > {min_specificity}? {is_specificity_valid}")
        
        if is_sensitivity_valid and is_specificity_valid:
            print(f"  Result: This threshold is valid.\n")
            valid_thresholds.append(threshold)
        else:
            print(f"  Result: This threshold is not valid.\n")

    # Find the maximum threshold among the valid ones
    if not valid_thresholds:
        print("No threshold found that satisfies both conditions.")
        # Exit gracefully if no answer is found
        sys.exit()

    max_valid_threshold = max(valid_thresholds)
    
    print("-------------------------------------------------")
    print(f"The valid thresholds are: {valid_thresholds}")
    print(f"The maximum Aortomesenteric diameter threshold that is >60% sensitive and >80% specific is: {max_valid_threshold} mm")
    print("-------------------------------------------------")

    return max_valid_threshold

# Execute the function and capture the final answer
final_answer = find_max_diameter_threshold()

# The final answer is enclosed in <<< >>>
print(f'<<<{final_answer}>>>')
