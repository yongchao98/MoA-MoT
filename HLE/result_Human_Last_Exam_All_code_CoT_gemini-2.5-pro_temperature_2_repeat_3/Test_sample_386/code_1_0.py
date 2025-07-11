def calculate_medium_composition():
    """
    Calculates the volumes of components for a cell culture medium based on
    the correct scenario from the multiple-choice question.

    The scenario involves a medium containing 10% FBS and 1% antibiotics.
    This script calculates the required volume of each component for a
    total volume of 500 mL.
    """

    total_volume_ml = 500
    fbs_percentage = 10
    antibiotic_percentage = 1

    # Calculate the volume for each component
    fbs_volume_ml = total_volume_ml * (fbs_percentage / 100)
    antibiotic_volume_ml = total_volume_ml * (antibiotic_percentage / 100)
    basal_medium_volume_ml = total_volume_ml - fbs_volume_ml - antibiotic_volume_ml

    # Print the explanation and the final equation
    print("To prepare 500 mL of fibroblast growth medium:")
    print(f"- Fetal Bovine Serum (FBS) at {fbs_percentage}%: {fbs_volume_ml:.0f} mL")
    print(f"- Antibiotics at {antibiotic_percentage}%: {antibiotic_volume_ml:.0f} mL")
    print(f"- Basal Medium: {basal_medium_volume_ml:.0f} mL")
    print("\nThe final medium composition equation is:")
    print(f"{total_volume_ml} mL Total = {basal_medium_volume_ml:.0f} mL Basal Medium + {fbs_volume_ml:.0f} mL ({fbs_percentage}%) FBS + {antibiotic_volume_ml:.0f} mL ({antibiotic_percentage}%) Antibiotics")

# Execute the function
calculate_medium_composition()