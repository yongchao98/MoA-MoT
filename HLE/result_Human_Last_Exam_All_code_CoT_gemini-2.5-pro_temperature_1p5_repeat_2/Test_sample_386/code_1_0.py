def prepare_culture_medium():
    """
    This function calculates the volumes of components needed to prepare
    cell culture medium based on the percentages given in the correct answer.
    It then prints the preparation recipe as an equation.
    """
    # Parameters from the correct option (C)
    fbs_percentage = 10  # %
    antibiotic_percentage = 1  # %

    # Standard total volume for a media bottle
    total_volume_ml = 500.0

    # Calculate the volume of each supplement
    fbs_volume_ml = total_volume_ml * (fbs_percentage / 100.0)
    antibiotic_volume_ml = total_volume_ml * (antibiotic_percentage / 100.0)

    # Calculate the volume of the base medium (e.g., DMEM)
    base_medium_volume_ml = total_volume_ml - fbs_volume_ml - antibiotic_volume_ml

    # Print the explanation and the final equation
    print("To prepare a 500 mL solution of fibroblast growth medium:")
    print(f"The recipe requires {fbs_percentage}% Fetal Bovine Serum (FBS) and {antibiotic_percentage}% Antibiotics.")
    print("\nCalculated Volumes:")
    print(f"Volume of Base Medium: {base_medium_volume_ml:.1f} mL")
    print(f"Volume of FBS ({fbs_percentage}%): {fbs_volume_ml:.1f} mL")
    print(f"Volume of Antibiotics ({antibiotic_percentage}%): {antibiotic_volume_ml:.1f} mL")
    print("\nFinal Preparation Equation:")
    # The following line explicitly outputs the numbers from the problem in the final equation.
    print(f"{base_medium_volume_ml:.1f} mL (Base Medium) + {fbs_volume_ml:.1f} mL ({fbs_percentage}% FBS) + {antibiotic_volume_ml:.1f} mL ({antibiotic_percentage}% Antibiotic) = {total_volume_ml:.1f} mL (Total Medium)")

# Execute the function
prepare_culture_medium()