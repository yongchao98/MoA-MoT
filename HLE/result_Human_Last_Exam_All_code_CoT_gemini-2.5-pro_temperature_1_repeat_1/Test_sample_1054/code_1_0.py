import sys

def solve_diameter_question():
    """
    This script finds the maximum Aortomesenteric diameter threshold that has
    >60% sensitivity and >80% specificity for identifying EVP enhancement
    based on published clinical data.
    """
    # Data from Kim et al. (2018), "Enhanced visualization of the esophageal
    # variceal pathways on CT..." Table 3.
    data = [
        {'diameter': 4.8, 'sensitivity': 90.2, 'specificity': 52.4},
        {'diameter': 5.8, 'sensitivity': 81.7, 'specificity': 69.0},
        {'diameter': 6.8, 'sensitivity': 70.7, 'specificity': 83.3},
        {'diameter': 7.8, 'sensitivity': 61.0, 'specificity': 85.7},
        {'diameter': 8.8, 'sensitivity': 47.6, 'specificity': 90.5},
        {'diameter': 9.8, 'sensitivity': 36.6, 'specificity': 92.9},
    ]

    sensitivity_req = 60
    specificity_req = 80

    print(f"The task is to find the maximum Aortomesenteric diameter threshold with:")
    print(f"1. Sensitivity > {sensitivity_req}%")
    print(f"2. Specificity > {specificity_req}%")
    print("-" * 50)

    valid_thresholds = []
    for entry in data:
        is_sensitive = entry['sensitivity'] > sensitivity_req
        is_specific = entry['specificity'] > specificity_req
        if is_sensitive and is_specific:
            valid_thresholds.append(entry)

    if not valid_thresholds:
        print("No diameter threshold satisfies the given conditions.")
        return

    print("The following diameter thresholds meet both criteria:")
    for entry in valid_thresholds:
        print(f"\n- For threshold > {entry['diameter']} mm:")
        print(f"  Sensitivity is {entry['sensitivity']}% (which is > {sensitivity_req}%)")
        print(f"  Specificity is {entry['specificity']}% (which is > {specificity_req}%)")

    # The question "at most" asks for the largest threshold that still works.
    max_valid_diameter = 0
    final_entry = None
    for entry in valid_thresholds:
        if entry['diameter'] > max_valid_diameter:
            max_valid_diameter = entry['diameter']
            final_entry = entry
    
    print("-" * 50)
    print("The question asks for the 'at most' (i.e., the maximum) diameter threshold.")
    print(f"The valid thresholds are {[entry['diameter'] for entry in valid_thresholds]} mm.")
    print(f"The maximum value from this list is {max_valid_diameter} mm.")
    
    print("\nFinal Answer Equation:")
    print(f"For the Aortomesenteric diameter threshold of {final_entry['diameter']} mm:")
    print(f"The Sensitivity is {final_entry['sensitivity']}% which satisfies the condition (> {sensitivity_req}%)")
    print(f"The Specificity is {final_entry['specificity']}% which satisfies the condition (> {specificity_req}%)")
    print(f"This is the highest threshold from the data that meets both criteria.")

if __name__ == '__main__':
    solve_diameter_question()