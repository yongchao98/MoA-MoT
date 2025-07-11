import math

def calculate_and_print_ferrite_level():
    """
    Calculates the approximate ferrite level for a given stainless steel composition
    and rounds it to the nearest 10.
    """
    # 1. Define the input equivalent values
    cr_eq = 39.0
    ni_eq = 29.0

    # Constants for the Schoefer formula
    cr_multiplier = 3.34
    ni_multiplier = 2.89
    subtractor = 27.57

    # 2. Calculate the approximate ferrite level using the formula
    ferrite_level = (cr_multiplier * cr_eq) - (ni_multiplier * ni_eq) - subtractor

    # 3. Round the result to the nearest 10
    # Python's round() rounds to the nearest even number for .5 cases (e.g., round(25, -1) is 20)
    # A more general method is to add 5, divide by 10, integer cast, then multiply by 10.
    # Or, even simpler: round(number / 10.0) * 10
    rounded_ferrite_level = int(round(ferrite_level / 10.0)) * 10

    # 4. Print the final equation with all numbers and the result
    print("Using the Schoefer formula to estimate the ferrite level:")
    print(f"Ferrite Level = ({cr_multiplier} * Chromium Eq.) - ({ni_multiplier} * Nickel Eq.) - {subtractor}")
    print("\nSubstituting the given values:")
    print(f"Ferrite Level = ({cr_multiplier} * {cr_eq}) - ({ni_multiplier} * {ni_eq}) - {subtractor}")

    # Calculate intermediate terms for the equation printout
    term1 = cr_multiplier * cr_eq
    term2 = ni_multiplier * ni_eq
    print(f"Ferrite Level = {term1:.2f} - {term2:.2f} - {subtractor}")
    print(f"Calculated Ferrite Level = {ferrite_level:.2f}")

    print(f"\nThe result rounded to the nearest 10 is: {rounded_ferrite_level}")

# Execute the function
calculate_and_print_ferrite_level()
<<<20>>>