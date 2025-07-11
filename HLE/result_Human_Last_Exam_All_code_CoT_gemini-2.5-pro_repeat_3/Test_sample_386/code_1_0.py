def prepare_culture_medium():
    """
    Calculates the volume of supplements needed to prepare a complete cell culture medium,
    based on the standard concentrations mentioned in the correct biological context.

    The context is the preparation of a medium for culturing corneal fibroblasts,
    which typically requires 10% Fetal Bovine Serum (FBS) and 1% antibiotics.
    """
    # Total volume of the complete medium to be prepared (in mL)
    total_volume_ml = 500

    # Percentages of supplements from the correct answer
    fbs_percentage = 10
    antibiotic_percentage = 1

    # Calculate the required volume for each supplement
    fbs_volume_ml = total_volume_ml * (fbs_percentage / 100.0)
    antibiotic_volume_ml = total_volume_ml * (antibiotic_percentage / 100.0)

    # Calculate the volume of the base medium required
    base_medium_ml = total_volume_ml - fbs_volume_ml - antibiotic_volume_ml

    print("To prepare {} mL of complete culture medium, you will need:\n".format(total_volume_ml))

    # Output the calculation for FBS
    print("Fetal Bovine Serum (FBS):")
    print("Equation: {} mL * {}% = {} mL".format(total_volume_ml, fbs_percentage, fbs_volume_ml))

    # Output the calculation for Antibiotics
    print("\nAntibiotics:")
    print("Equation: {} mL * {}% = {} mL".format(total_volume_ml, antibiotic_percentage, antibiotic_volume_ml))
    
    print("\nBase Medium: {} mL".format(base_medium_ml))

# Run the calculation
prepare_culture_medium()