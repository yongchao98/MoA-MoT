import math

def calculate_and_display_ferrite_level():
    """
    Calculates the approximate ferrite level for a stainless steel given
    its chromium and nickel equivalents.
    """
    # Given equivalent values
    ni_eq = 29.0
    cr_eq = 39.0

    # Coefficients for a simplified linear approximation formula:
    # Ferrite % = C1 * Cr_eq - C2 * Ni_eq - C3
    coeff_cr = 4.0
    coeff_ni = 2.5
    constant = 30.0

    # Calculate the approximate ferrite percentage
    ferrite_percent = coeff_cr * cr_eq - coeff_ni * ni_eq - constant

    # As requested, output the full equation with the numbers
    print("Using the approximate formula: Ferrite % = (4.0 * Cr_eq) - (2.5 * Ni_eq) - 30")
    print(f"Calculation: Ferrite % = ({coeff_cr} * {cr_eq}) - ({coeff_ni} * {ni_eq}) - {constant}")
    print(f"Result before rounding: {ferrite_percent:.1f}")

    # Round the result to the nearest 10
    # For a value like 53.5, it is 3.5 away from 50 and 6.5 away from 60.
    # Therefore, the nearest 10 is 50.
    # The `round(number / 10) * 10` logic correctly handles this.
    rounded_ferrite = round(ferrite_percent / 10) * 10
    
    print(f"\nThe calculated ferrite level rounded to the nearest 10 is: {rounded_ferrite}")

# Run the calculation
calculate_and_display_ferrite_level()