import math

def calculate_ferrite_level():
    """
    Calculates the approximate ferrite level for a stainless steel using Hull's formula
    and rounds it to the nearest 10.
    """
    # Given values
    ni_eq = 29.0
    cr_eq = 39.0
    
    # Hull's formula for Ferrite Number (FN) which approximates ferrite percentage
    # FN = 3.34 * Cr_eq - 2.46 * Ni_eq - 28.6
    
    term1 = 3.34 * cr_eq
    term2 = 2.46 * ni_eq
    constant = 28.6
    
    ferrite_level = term1 - term2 - constant
    
    # Print the equation with the substituted values
    print(f"Calculation based on Hull's formula:")
    print(f"Ferrite Level = 3.34 * {cr_eq} - 2.46 * {ni_eq} - {constant}")
    print(f"Calculated Value: {ferrite_level:.2f}")

    # Round the result to the nearest 10
    # The formula is to divide by 10, round to the nearest integer, then multiply by 10.
    rounded_ferrite_level = int(round(ferrite_level / 10.0) * 10)
    
    print(f"Rounded to the nearest 10: {rounded_ferrite_level}")
    
    # Final answer format
    print(f"<<<{rounded_ferrite_level}>>>")

calculate_ferrite_level()