import math

def calculate_ferrite_level():
    """
    Calculates the approximate ferrite level for a stainless steel
    using an approximation of the WRC-1992 diagram.
    """
    # 1. Define the input values
    cr_eq = 39.0
    ni_eq = 29.0

    # 2. Use Kotecki's formula to calculate the Ferrite Number (FN)
    # Breaking down the formula for clarity and to show each number
    c1 = -140.2
    c2 = 8.23
    c3 = -4.93
    c4 = 0.117
    c5 = 0.038
    c6 = -0.0049

    term1 = c1
    term2 = c2 * cr_eq
    term3 = c3 * ni_eq
    term4 = c4 * (cr_eq**2)
    term5 = c5 * (ni_eq**2)
    term6 = c6 * (cr_eq**2) * ni_eq

    ferrite_number = term1 + term2 + term3 + term4 + term5 + term6

    # 4. Round the result to the nearest 10
    rounded_ferrite = round(ferrite_number / 10) * 10
    
    # 5. Print the full equation with the numbers used
    print("The calculation is based on an approximation of the WRC-1992 diagram.")
    print(f"The formula used is: FN = {c1} + ({c2} * Cr_eq) + ({c3} * Ni_eq) + ({c4} * Cr_eq^2) + ({c5} * Ni_eq^2) + ({c6} * Cr_eq^2 * Ni_eq)")
    print("\nPlugging in the values:")
    print(f"FN = {c1} + ({c2} * {cr_eq}) + ({c3} * {ni_eq}) + ({c4} * {cr_eq}^2) + ({c5} * {ni_eq}^2) + ({c6} * {cr_eq}^2 * {ni_eq})")
    
    # 6. Print the result
    print(f"\nCalculated Ferrite Number (FN): {ferrite_number:.2f}")
    print(f"Approximate ferrite level rounded to the nearest 10 is: {rounded_ferrite}")

if __name__ == '__main__':
    calculate_ferrite_level()