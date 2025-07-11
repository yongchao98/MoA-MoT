def find_different_imprints():
    """
    Identifies and prints the different numbers on Teva brand fluoxetine
    10mg and 20mg capsules.
    """
    # Numerical imprints on the capsules, excluding non-numeric parts
    imprint_10mg = "3105"
    imprint_20mg = "3106"

    # Create sets of the individual numbers (digits) for easy comparison
    set_10mg = set(imprint_10mg)
    set_20mg = set(imprint_20mg)

    # Find the numbers that are not present in both sets (the symmetric difference)
    different_numbers = sorted(list(set_10mg.symmetric_difference(set_20mg)))

    # Output the results as requested
    print("The numbers on the 10mg capsule are: " + ", ".join(list(imprint_10mg)))
    print("The numbers on the 20mg capsule are: " + ", ".join(list(imprint_20mg)))
    print("\nThe different numbers imprinted on the capsules are: " + ", ".join(different_numbers))

find_different_imprints()