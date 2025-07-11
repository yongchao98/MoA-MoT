import re

def find_imprint_differences():
    """
    Identifies and prints the different numbers on Teva fluoxetine 10mg and 20mg capsules.
    """
    # The full imprints are typically "TEVA 3108" and "TEVA 3109".
    # We will focus on the numerical parts as requested.
    imprint_10mg_num = "3108"
    imprint_20mg_num = "3109"

    # Print the individual numbers for each capsule's imprint
    print(f"The Teva fluoxetine 10mg capsule has the numerical imprint: {' '.join(imprint_10mg_num)}")
    print(f"The Teva fluoxetine 20mg capsule has the numerical imprint: {' '.join(imprint_20mg_num)}")

    # Use sets to find the digits that are different
    set_10mg = set(imprint_10mg_num)
    set_20mg = set(imprint_20mg_num)

    # The symmetric difference will give us the elements present in one set but not the other.
    different_numbers = sorted(list(set_10mg.symmetric_difference(set_20mg)))

    print(f"\nComparing the two numerical imprints, the numbers that are different are {different_numbers[0]} and {different_numbers[1]}.")

find_imprint_differences()