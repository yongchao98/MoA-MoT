import math

def calculate_ferrite_level():
    """
    Calculates the approximate ferrite level for a stainless steel
    using Hull's formula and rounds it to the nearest 10.
    """
    # Given equivalent values
    cr_eq = 39
    ni_eq = 29

    # Hull's formula for approximate ferrite percentage
    # Ferrite % = -40.66 + 4.47 * Creq - 3.44 * Nieq
    ferrite_percent = -40.66 + 4.47 * cr_eq - 3.44 * ni_eq

    # Round the result to the nearest 10
    # The round(number, -1) function rounds to the nearest 10
    rounded_ferrite = int(round(ferrite_percent, -1))

    # Print the equation and the final result
    print(f"Calculating ferrite percentage using Hull's formula:")
    print(f"Ferrite % = -40.66 + (4.47 * {cr_eq}) - (3.44 * {ni_eq})")
    print(f"Ferrite % = -40.66 + {4.47 * cr_eq} - {3.44 * ni_eq}")
    print(f"Ferrite % = {ferrite_percent:.2f}")
    print(f"\nRounding {ferrite_percent:.2f} to the nearest 10 results in {rounded_ferrite}.")
    print(f"The approximate ferrite level is {rounded_ferrite}.")

calculate_ferrite_level()
<<<30>>>