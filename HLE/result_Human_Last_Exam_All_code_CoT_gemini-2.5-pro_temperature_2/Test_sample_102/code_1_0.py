import math

def calculate_ferrite_level():
    """
    Calculates the approximate ferrite level for a stainless steel using the
    Kotecki formula based on the WRC-1992 diagram and rounds it.
    """
    # Input values from the user's request
    cr_eq = 39.0
    ni_eq = 29.0

    # Constants from Kotecki's formula
    c1 = 5.26
    c2 = 0.74
    c3 = 9.49

    # Calculate the ferrite level
    # The result is the Ferrite Number (FN), which is approximately the ferrite percentage
    ferrite_level = c1 * (cr_eq - (c2 * ni_eq) - c3)

    # Round the calculated ferrite level to the nearest 10
    rounded_ferrite_level = round(ferrite_level / 10) * 10

    # Print the equation with all the numbers and the unrounded result
    # We show each number that goes into the final equation
    print(f"The ferrite level is calculated using the equation: {c1} * ({cr_eq} - {c2} * {ni_eq} - {c3})")
    print(f"Calculation result (unrounded): {ferrite_level:.1f}%")
    print("-" * 30)

    # Print the final, rounded answer
    print(f"The approximate ferrite level, rounded to the nearest 10, is: {rounded_ferrite_level}")

if __name__ == "__main__":
    calculate_ferrite_level()