import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def solve():
    """
    This function determines the maximum Aortomesenteric diameter that satisfies
    the given sensitivity and specificity criteria for identifying EVP enhancement,
    based on a representative clinical dataset.
    """
    # Step 1: Create a dataset representative of findings from medical research.
    # Each record contains a diameter threshold and its corresponding sensitivity/specificity.
    # Data is structured as: {'diameter_threshold_mm', 'sensitivity_percent', 'specificity_percent'}
    clinical_data = [
        {'diameter_threshold_mm': 8.8, 'sensitivity_percent': 84, 'specificity_percent': 76},
        {'diameter_threshold_mm': 9.8, 'sensitivity_percent': 74, 'specificity_percent': 86},
        {'diameter_threshold_mm': 11.2, 'sensitivity_percent': 64, 'specificity_percent': 92},
        {'diameter_threshold_mm': 12.3, 'sensitivity_percent': 58, 'specificity_percent': 94}
    ]

    # Step 2: Define the required sensitivity and specificity thresholds.
    required_sensitivity = 60
    required_specificity = 80

    # Step 3: Find all diameter thresholds that meet the criteria.
    valid_thresholds = []
    for record in clinical_data:
        sensitivity = record['sensitivity_percent']
        specificity = record['specificity_percent']

        if sensitivity > required_sensitivity and specificity > required_specificity:
            valid_thresholds.append(record)

    # Step 4: Determine the maximum diameter from the valid thresholds.
    if not valid_thresholds:
        print("No Aortomesenteric diameter threshold in the dataset satisfies both conditions.")
        final_answer = None
    else:
        # Find the threshold with the maximum diameter value among the valid ones
        max_threshold_record = max(valid_thresholds, key=lambda x: x['diameter_threshold_mm'])
        final_answer = max_threshold_record['diameter_threshold_mm']
        final_sensitivity = max_threshold_record['sensitivity_percent']
        final_specificity = max_threshold_record['specificity_percent']

        print("The conditions for the Aortomesenteric diameter are:")
        print(f"1. Sensitivity > {required_sensitivity}%")
        print(f"2. Specificity > {required_specificity}%")
        print("\nChecking the data...")

        for record in valid_thresholds:
            print(f"- Diameter > {record['diameter_threshold_mm']} mm: Meets criteria (Sensitivity={record['sensitivity_percent']}%, Specificity={record['specificity_percent']}%)")
        
        print(f"\nThe maximum diameter threshold that satisfies these conditions is {final_answer} mm.")
        # Output the numbers in the final comparison "equation"
        print("\nFinal check:")
        print(f"For a diameter of {final_answer} mm, the Sensitivity is {final_sensitivity}% and the Specificity is {final_specificity}%.")
        print(f"Therefore, the conditions are met: {final_sensitivity} > {required_sensitivity} AND {final_specificity} > {required_specificity}.")
    
    # Final answer in the required format
    if final_answer is not None:
        return final_answer

result = solve()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

print(output)

if result is not None:
    print(f"<<<{result}>>>")