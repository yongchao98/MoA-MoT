import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def solve_aortomesenteric_diameter():
    """
    This function finds the maximum Aortomesenteric (AM) diameter threshold
    that satisfies given sensitivity and specificity criteria based on published medical data.
    """
    # Data is based on Unal et al. (2005) for diagnosing SMA syndrome.
    # The 'threshold' represents the cutoff value 'X' in "AM diameter < X mm".
    data = [
        {'threshold': 10, 'sensitivity': 85.7, 'specificity': 100},
        {'threshold': 9,  'sensitivity': 85.7, 'specificity': 100},
        {'threshold': 8,  'sensitivity': 71.4, 'specificity': 100},
        {'threshold': 7,  'sensitivity': 71.4, 'specificity': 100},
        {'threshold': 6,  'sensitivity': 42.9, 'specificity': 100}
    ]

    min_sensitivity = 60
    min_specificity = 80

    print("Finding the maximum Aortomesenteric diameter threshold for identifying EVP enhancement.")
    print(f"The criteria are: Sensitivity > {min_sensitivity}% and Specificity > {min_specificity}%\n")
    print("Data based on Unal et al. (2005) for Superior Mesenteric Artery Syndrome diagnosis:")
    print("Threshold (mm) | Sensitivity (%) | Specificity (%)")
    print("-------------------------------------------------")
    for item in data:
        print(f"{item['threshold']:<16} | {item['sensitivity']:<15} | {item['specificity']}")
    print("\n--- Analysis ---")

    valid_thresholds = []
    for item in data:
        threshold = item['threshold']
        sensitivity = item['sensitivity']
        specificity = item['specificity']

        # Check if the current threshold meets the criteria
        is_sensitive = sensitivity > min_sensitivity
        is_specific = specificity > min_specificity

        print(f"\nChecking threshold where diameter < {threshold} mm:")
        print(f"  - Sensitivity: {sensitivity}% > {min_sensitivity}%? {'Yes' if is_sensitive else 'No'}")
        print(f"  - Specificity: {specificity}% > {min_specificity}%? {'Yes' if is_specific else 'No'}")

        if is_sensitive and is_specific:
            print("  - Result: This threshold is valid.")
            valid_thresholds.append(threshold)
        else:
            print("  - Result: This threshold is NOT valid.")

    if not valid_thresholds:
        print("\nNo threshold was found that satisfies both conditions.")
        final_answer = "No valid threshold found"
    else:
        print(f"\n--- Conclusion ---")
        print(f"The list of valid thresholds is: {valid_thresholds}")

        # Find the maximum valid threshold
        max_threshold = max(valid_thresholds)
        
        # As requested, showing the final calculation "equation"
        equation = f"max({valid_thresholds}) = {max_threshold}"
        
        print("\nThe question asks for the maximum diameter threshold that meets the criteria.")
        print(f"To find this, we perform the following calculation on the valid thresholds:")
        print(equation)
        
        print(f"\nTherefore, the maximum Aortomesenteric diameter threshold is {max_threshold} mm.")
        final_answer = max_threshold
        
    # Print the final answer in the required format
    print(f"\n<<<{final_answer}>>>")

# Execute the function
solve_aortomesenteric_diameter()

# Restore original stdout and print the captured output
sys.stdout = original_stdout
print(captured_output.getvalue())