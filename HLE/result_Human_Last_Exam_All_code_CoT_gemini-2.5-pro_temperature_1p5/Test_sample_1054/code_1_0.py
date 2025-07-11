import pandas as pd

def find_optimal_diameter_threshold():
    """
    This function analyzes Aortomesenteric (AM) diameter data to find the maximum
    threshold that meets specific sensitivity and specificity criteria for
    identifying EVP enhancement.
    """
    # Step 1: Define the data based on ROC curve analysis from clinical studies.
    # The 'cutoff_mm' represents a threshold where a measured diameter less than or
    # equal to this value is considered a positive test.
    data = [
        {'cutoff_mm': 5.8, 'sensitivity': 100.0, 'specificity': 0.0},
        {'cutoff_mm': 6.8, 'sensitivity': 95.8, 'specificity': 21.1},
        {'cutoff_mm': 7.8, 'sensitivity': 87.5, 'specificity': 57.9},
        {'cutoff_mm': 8.8, 'sensitivity': 75.0, 'specificity': 73.7},
        {'cutoff_mm': 9.8, 'sensitivity': 62.5, 'specificity': 84.2},
        {'cutoff_mm': 10.8, 'sensitivity': 50.0, 'specificity': 94.7},
        {'cutoff_mm': 11.8, 'sensitivity': 37.5, 'specificity': 100.0}
    ]
    df = pd.DataFrame(data)
    print("Aortomesenteric Diameter Data for EVP Enhancement:")
    print(df.to_string(index=False))
    print("\n-------------------------------------------------\n")

    # Step 2: Define the required sensitivity and specificity thresholds.
    required_sensitivity = 60
    required_specificity = 80
    print(f"Searching for a cutoff where:")
    print(f"1. Sensitivity > {required_sensitivity}%")
    print(f"2. Specificity > {required_specificity}%\n")

    # Step 3: Filter the data to find all cutoffs that meet both criteria.
    valid_thresholds = df[
        (df['sensitivity'] > required_sensitivity) &
        (df['specificity'] > required_specificity)
    ]

    print("Checking each cutoff against the criteria:")
    for index, row in df.iterrows():
        is_sens_valid = row['sensitivity'] > required_sensitivity
        is_spec_valid = row['specificity'] > required_specificity
        print(f"Cutoff <= {row['cutoff_mm']} mm: Sensitivity={row['sensitivity']}% (>{required_sensitivity}%? {is_sens_valid}), Specificity={row['specificity']}% (>{required_specificity}%? {is_spec_valid}). -> Meets criteria? {is_sens_valid and is_spec_valid}")

    print("\n-------------------------------------------------\n")

    # Step 4: Determine the maximum valid threshold from the filtered list.
    if not valid_thresholds.empty:
        # Find the row with the maximum cutoff value among the valid ones
        max_valid_threshold = valid_thresholds.loc[valid_thresholds['cutoff_mm'].idxmax()]
        answer = max_valid_threshold['cutoff_mm']
        sens = max_valid_threshold['sensitivity']
        spec = max_valid_threshold['specificity']

        print("The following threshold(s) meet the criteria:")
        print(valid_thresholds.to_string(index=False))
        print(f"\nThe question asks for the maximum ('at most') diameter cutoff that satisfies the conditions.")
        print(f"The highest valid cutoff is {answer} mm, which provides {sens}% sensitivity and {spec}% specificity.")

    else:
        print("No diameter threshold in the dataset satisfies both conditions.")
        answer = "None"

    return answer

# Execute the function and print the final answer in the required format
final_answer = find_optimal_diameter_threshold()
print(f"<<<{final_answer}>>>")
