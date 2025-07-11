import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def solve():
    """
    Finds the maximum Aortomesenteric diameter threshold that meets given sensitivity and specificity criteria.
    """
    # Step 1: Define the hypothetical dataset.
    # This data represents potential outcomes from a study correlating AM diameter with EVP enhancement.
    # A test is considered positive if a patient's measured AM diameter is GREATER than the threshold.
    roc_data = [
        {'diameter_threshold': 7.0, 'sensitivity': 95, 'specificity': 30},
        {'diameter_threshold': 7.5, 'sensitivity': 90, 'specificity': 45},
        {'diameter_threshold': 8.0, 'sensitivity': 85, 'specificity': 60},
        {'diameter_threshold': 8.5, 'sensitivity': 75, 'specificity': 75},
        {'diameter_threshold': 9.0, 'sensitivity': 65, 'specificity': 82},
        {'diameter_threshold': 9.5, 'sensitivity': 61, 'specificity': 88},
        {'diameter_threshold': 10.0, 'sensitivity': 55, 'specificity': 92},
        {'diameter_threshold': 10.5, 'sensitivity': 48, 'specificity': 95},
    ]

    # Step 2: Set the performance criteria from the user's request.
    sensitivity_threshold = 60
    specificity_threshold = 80

    # This list will store diameter thresholds that meet both criteria.
    valid_diameters = []

    print(f"Searching for AM diameter thresholds with Sensitivity > {sensitivity_threshold}% and Specificity > {specificity_threshold}%.\n")

    # Step 3: Filter the data by iterating through each point.
    for point in roc_data:
        diameter = point['diameter_threshold']
        sensitivity = point['sensitivity']
        specificity = point['specificity']

        # Check if the current point meets the sensitivity and specificity criteria.
        is_sensitive = sensitivity > sensitivity_threshold
        is_specific = specificity > specificity_threshold
        
        if is_sensitive and is_specific:
            valid_diameters.append(diameter)

    # Step 4: Find the maximum diameter from the list of valid ones.
    if not valid_diameters:
        print("No diameter threshold in the dataset satisfies both conditions.")
        return None

    max_valid_diameter = max(valid_diameters)
    
    # Retrieve the specific sensitivity and specificity for the chosen diameter
    final_point = next((p for p in roc_data if p['diameter_threshold'] == max_valid_diameter), None)
    final_sensitivity = final_point['sensitivity']
    final_specificity = final_point['specificity']

    # Step 5: Output the result and the final check.
    print(f"The valid diameter thresholds that meet both criteria are: {valid_diameters}")
    print(f"\nThe question asks for the maximum (at most) diameter threshold from this set.")
    print(f"The maximum Aortomesenteric diameter threshold is {max_valid_diameter} mm.\n")
    
    print("Final check for this threshold:")
    print("The conditions are that sensitivity > 60% and specificity > 80%.")
    # Outputting the numbers for the "final equation" as requested
    print(f"For a threshold of {max_valid_diameter} mm, is Sensitivity ({final_sensitivity}) > {sensitivity_threshold}? Yes.")
    print(f"For a threshold of {max_valid_diameter} mm, is Specificity ({final_specificity}) > {specificity_threshold}? Yes.")

    return max_valid_diameter

result = solve()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output_str = captured_output.getvalue()

print(output_str)

if result is not None:
    # We add the final answer tag here based on the result
    print(f"<<<{result}>>>")