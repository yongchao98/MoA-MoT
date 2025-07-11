import math

def calculate_ferrite_level():
    """
    Calculates the approximate ferrite level for a stainless steel
    based on its nickel and chromium equivalents.
    """
    # Given values
    nickel_equivalent = 29
    chromium_equivalent = 39

    # An empirical formula to approximate the Ferrite Number (FN)
    # FN is considered roughly equal to the ferrite percentage.
    # Formula: FN = 3.34 * Cr_eq - 2.46 * Ni_eq - 28.6
    ferrite_level = 3.34 * chromium_equivalent - 2.46 * nickel_equivalent - 28.6

    # Round the calculated ferrite level to the nearest 10
    rounded_ferrite_level = int(round(ferrite_level / 10.0) * 10)

    # Print the equation with all the numbers and the final rounded answer
    print(f"3.34 * {chromium_equivalent} - 2.46 * {nickel_equivalent} - 28.6 = {rounded_ferrite_level} (Result {ferrite_level:.2f} rounded to the nearest 10)")

calculate_ferrite_level()