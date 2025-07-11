import math

def calculate_ferrite_level():
    """
    Calculates the approximate ferrite level for a stainless steel based on its
    nickel and chromium equivalents, then rounds the result to the nearest 10.
    """
    # Given equivalent values
    ni_eq = 29.0
    cr_eq = 39.0

    # Constants for a common empirical formula to approximate ferrite level
    c1 = 3.34
    c2 = 2.46
    c3 = 28.6

    print(f"Calculating the approximate ferrite level for a stainless steel with:")
    print(f"Nickel Equivalent (Ni_eq) = {ni_eq}")
    print(f"Chromium Equivalent (Cr_eq) = {cr_eq}\n")

    print(f"Using the formula: Ferrite Level = ({c1} * Cr_eq) - ({c2} * Ni_eq) - {c3}")

    # Calculate the ferrite level
    ferrite_level = (c1 * cr_eq) - (c2 * ni_eq) - c3

    # Ferrite level cannot be negative
    if ferrite_level < 0:
        ferrite_level = 0.0

    # Round the result to the nearest 10
    rounded_ferrite_level = int(round(ferrite_level, -1))

    print("\nCalculation steps:")
    # The final equation with all numbers plugged in
    print(f"Ferrite Level = ({c1} * {cr_eq}) - ({c2} * {ni_eq}) - {c3}")
    print(f"Ferrite Level = {c1 * cr_eq:.2f} - {c2 * ni_eq:.2f} - {c3}")
    print(f"Calculated Ferrite Level = {ferrite_level:.2f}\n")
    print(f"Rounding to the nearest 10, the final approximate ferrite level is: {rounded_ferrite_level}")


if __name__ == '__main__':
    calculate_ferrite_level()