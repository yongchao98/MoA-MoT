def find_different_imprints():
    """
    Identifies and prints the different numbers on Teva fluoxetine capsules.
    """
    # Numerical imprints for Teva fluoxetine 10mg and 20mg capsules
    imprint_10mg = "3108"
    imprint_20mg = "3109"

    # Convert the imprints to sets of characters to find the difference
    set_10mg = set(imprint_10mg)
    set_20mg = set(imprint_20mg)

    # The symmetric difference finds elements that are in one set, but not in both
    different_numbers = sorted(list(set_10mg.symmetric_difference(set_20mg)))

    # Print the findings in a descriptive sentence
    print(f"The number imprinted on the 10mg capsule is {imprint_10mg}.")
    print(f"The number imprinted on the 20mg capsule is {imprint_20mg}.")
    print(f"The different numbers between them are {different_numbers[0]} and {different_numbers[1]}.")

find_different_imprints()