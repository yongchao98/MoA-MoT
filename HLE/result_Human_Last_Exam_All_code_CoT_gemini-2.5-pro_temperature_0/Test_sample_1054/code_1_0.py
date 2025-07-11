import pandas as pd
import io

def solve():
    """
    This function finds the maximum Aortomesenteric diameter that satisfies
    both >60% sensitivity and >80% specificity for identifying EVP enhancement,
    based on a given dataset.
    """
    # Step 1: Establish a dataset.
    # In a real-world application, this data would be loaded from a file (e.g., CSV)
    # from a clinical study. Here, it's defined directly in the code for demonstration.
    data = """
    AM_Diameter_Cutoff_mm,Sensitivity_Percent,Specificity_Percent
    5.0,95,40
    6.0,90,55
    7.0,85,70
    8.0,75,82
    9.0,65,88
    10.0,55,92
    11.0,45,95
    12.0,30,98
    """
    df = pd.read_csv(io.StringIO(data))

    # Step 2: Define the criteria thresholds.
    sensitivity_threshold = 60
    specificity_threshold = 80

    # Initialize a variable to store the maximum valid diameter.
    max_valid_diameter = None
    final_sensitivity = None
    final_specificity = None

    print("Evaluating each Aortomesenteric diameter cutoff:")
    # Step 3: Iterate through the dataset and evaluate each cutoff.
    for index, row in df.iterrows():
        diameter = row['AM_Diameter_Cutoff_mm']
        sensitivity = row['Sensitivity_Percent']
        specificity = row['Specificity_Percent']

        # Check if the current cutoff meets the criteria.
        is_sensitive_enough = sensitivity > sensitivity_threshold
        is_specific_enough = specificity > specificity_threshold

        print(f"- Checking Diameter {diameter} mm: Sensitivity={sensitivity}%, Specificity={specificity}%.")

        if is_sensitive_enough and is_specific_enough:
            print(f"  -> Condition met: Sensitivity ({sensitivity}%) > {sensitivity_threshold}% AND Specificity ({specificity}%) > {specificity_threshold}%.")
            # Step 4: If criteria are met, update the maximum valid diameter.
            if max_valid_diameter is None or diameter > max_valid_diameter:
                max_valid_diameter = diameter
                final_sensitivity = sensitivity
                final_specificity = specificity
        else:
            print(f"  -> Condition not met.")


    # Step 5: Output the final result.
    print("\n" + "="*50)
    if max_valid_diameter is not None:
        print("Final Answer:")
        # The "final equation" is showing the numbers that satisfy the condition for the maximum diameter found.
        print(f"The maximum Aortomesenteric diameter that is >{sensitivity_threshold}% sensitive and >{specificity_threshold}% specific is {max_valid_diameter} mm.")
        print(f"This is based on the values: Sensitivity = {final_sensitivity}% and Specificity = {final_specificity}%.")
    else:
        print("No Aortomesenteric diameter in the dataset satisfies both conditions.")
    print("="*50)

solve()
<<<9.0>>>