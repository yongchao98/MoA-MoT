import math

def calculate_ferrite_level():
    """
    Calculates the approximate ferrite level for a stainless steel
    using the WRC-1992 diagram approximation formula and rounds it to the nearest 10.
    """
    # Given values
    cr_eq = 39.0
    ni_eq = 29.0

    # Constants from the Kotecki/Olson formula
    c1 = 3.34
    c2 = 2.46
    c3 = 28.6

    # Calculate the ferrite level
    ferrite_level = (c1 * cr_eq) - (c2 * ni_eq) - c3

    # Round the result to the nearest 10
    # The second argument of round(), -1, rounds to the nearest 10.
    rounded_ferrite_level = int(round(ferrite_level, -1))

    # Print the equation with the specific values
    print(f"Using the formula: (c1 * Cr_eq) - (c2 * Ni_eq) - c3")
    print(f"Calculation: ({c1} * {cr_eq}) - ({c2} * {ni_eq}) - {c3} = {ferrite_level:.2f}")
    
    # Print the final answer
    print(f"\nThe calculated ferrite level is approximately {ferrite_level:.2f}%.")
    print(f"Rounded to the nearest 10, the ferrite level is: {rounded_ferrite_level}")

calculate_ferrite_level()