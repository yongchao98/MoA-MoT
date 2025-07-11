import math

def calculate_and_print_ferrite_level():
    """
    Calculates the approximate ferrite level for a stainless steel
    using the Hull formula and rounds the result to the nearest 10.
    """
    # Given equivalents
    cr_eq = 39
    ni_eq = 29

    # Constants for the Hull formula
    c1 = 3.34
    c2 = 2.46
    c3 = 28.6

    # Calculate the intermediate terms
    term1 = c1 * cr_eq
    term2 = c2 * ni_eq

    # Calculate the unrounded ferrite level
    ferrite_level = term1 - term2 - c3

    # Round the final result to the nearest 10
    # To round to the nearest 10, we divide by 10, round to the nearest integer, then multiply by 10.
    rounded_ferrite_level = round(ferrite_level / 10) * 10

    # Print the equation and the result
    print("Using the Hull formula to estimate ferrite level:")
    print(f"Ferrite Level = ({c1} * Cr_eq) - ({c2} * Ni_eq) - {c3}")
    print("\nSubstituting the given values:")
    print(f"Ferrite Level = ({c1} * {cr_eq}) - ({c2} * {ni_eq}) - {c3}")
    print(f"Ferrite Level = {term1:.2f} - {term2:.2f} - {c3}")
    print(f"\nCalculated Ferrite Level (unrounded): {ferrite_level:.2f}")
    print(f"Approximate Ferrite Level (rounded to the nearest 10): {rounded_ferrite_level}")

# Execute the function
calculate_and_print_ferrite_level()