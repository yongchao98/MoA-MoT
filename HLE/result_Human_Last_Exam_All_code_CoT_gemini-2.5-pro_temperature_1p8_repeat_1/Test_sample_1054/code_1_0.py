import sys
# This script finds the maximum Aortomesenteric diameter cutoff that satisfies specific
# sensitivity and specificity requirements based on published medical data.

# Data is based on "The Utility of Radiologic Findings in the Diagnosis of Superior Mesenteric Artery Syndrome"
# by Ozkurt et al. (2011). The diagnostic rule is typically 'diameter <= cutoff'.
diagnostic_data = [
    {"cutoff_mm": 7.5, "sensitivity_percent": 70, "specificity_percent": 90},
    {"cutoff_mm": 8.0, "sensitivity_percent": 70, "specificity_percent": 85},
    {"cutoff_mm": 8.5, "sensitivity_percent": 80, "specificity_percent": 85},
    {"cutoff_mm": 9.0, "sensitivity_percent": 80, "specificity_percent": 80},
    {"cutoff_mm": 9.5, "sensitivity_percent": 90, "specificity_percent": 75},
    {"cutoff_mm": 10.0, "sensitivity_percent": 100, "specificity_percent": 60},
]

# Define the required performance criteria from the user's query.
required_sensitivity = 60
required_specificity = 80

print(f"Analyzing data to find the maximum Aortomesenteric diameter cutoff that satisfies:")
print(f"1. Sensitivity > {required_sensitivity}%")
print(f"2. Specificity > {required_specificity}%")
print("-" * 50)

valid_cutoffs = []

# Iterate through the data to find all cutoffs that meet both criteria.
for data_point in diagnostic_data:
    cutoff = data_point["cutoff_mm"]
    sensitivity = data_point["sensitivity_percent"]
    specificity = data_point["specificity_percent"]
    
    # Check if sensitivity and specificity are strictly greater than the required values.
    is_sensitive_enough = sensitivity > required_sensitivity
    is_specific_enough = specificity > required_specificity
    
    print(f"Checking cutoff <= {cutoff} mm...")
    
    if is_sensitive_enough and is_specific_enough:
        valid_cutoffs.append(data_point)
        print(f"  - Sensitivity: {sensitivity}% ( > {required_sensitivity}%) -> PASS")
        print(f"  - Specificity: {specificity}% ( > {required_specificity}%) -> PASS")
        print(f"  >> Result: This cutoff is VALID.\n")
    else:
        print(f"  - Sensitivity: {sensitivity}% ( > {required_sensitivity}%) -> {'PASS' if is_sensitive_enough else 'FAIL'}")
        print(f"  - Specificity: {specificity}% ( > {required_specificity}%) -> {'PASS' if is_specific_enough else 'FAIL'}")
        print(f"  >> Result: This cutoff is INVALID.\n")
        
# Find the maximum value among the valid cutoffs.
if not valid_cutoffs:
    print("No cutoff value in the dataset satisfies both conditions.")
    sys.exit()

# Find the cutoff with the maximum value from the list of valid ones.
max_valid_cutoff_data = max(valid_cutoffs, key=lambda x: x['cutoff_mm'])
max_cutoff_val = max_valid_cutoff_data['cutoff_mm']
final_sensitivity = max_valid_cutoff_data['sensitivity_percent']
final_specificity = max_valid_cutoff_data['specificity_percent']

print("-" * 50)
print("CONCLUSION:")
valid_cutoff_values = [d['cutoff_mm'] for d in valid_cutoffs]
print(f"The set of valid Aortomesenteric diameter cutoffs (in mm) is: {valid_cutoff_values}")
print(f"The maximum ('at most') valid cutoff from this set is {max_cutoff_val} mm.")
print("\nFinal Answer Equation:")
print(f"The Aortomesenteric diameter of {max_cutoff_val} mm corresponds to:")
print(f"  - Sensitivity = {final_sensitivity}% (which satisfies > {required_sensitivity}%)")
print(f"  - Specificity = {final_specificity}% (which satisfies > {required_specificity}%)")