import math

def calculate_ferrite_level():
    """
    Calculates the approximate ferrite level for a stainless steel
    based on its nickel and chromium equivalents.
    """
    # Step 1: Define the given equivalent values
    ni_eq = 29.0
    cr_eq = 39.0
    
    # Coefficients from the WRC-1992 simplified formula
    coeff_cr = 3.34
    coeff_ni = 2.46
    constant = -28.6

    print("To approximate the ferrite level, a formula derived from the WRC-1992 diagram is used.")
    print(f"The formula is: Ferrite Level = {coeff_cr} * Cr_eq - {coeff_ni} * Ni_eq + ({constant})")
    print("\nGiven values:")
    print(f"Nickel Equivalent (Ni_eq) = {ni_eq}")
    print(f"Chromium Equivalent (Cr_eq) = {cr_eq}")
    
    # Step 2: Calculate the ferrite level using the formula
    ferrite_level = coeff_cr * cr_eq + coeff_ni * -ni_eq + constant
    
    # Step 3: Print the calculation steps
    print("\nCalculation:")
    print(f"Ferrite Level = {coeff_cr} * {cr_eq} - {coeff_ni} * {ni_eq} - {abs(constant)}")
    print(f"Ferrite Level = {coeff_cr * cr_eq:.2f} - {coeff_ni * ni_eq:.2f} - {abs(constant)}")
    print(f"Ferrite Level = {ferrite_level:.2f}")

    # Step 4: Round the result to the nearest 10
    rounded_ferrite_level = int(round(ferrite_level / 10.0)) * 10
    
    print(f"\nThe calculated ferrite level is approximately {ferrite_level:.2f}%.")
    print(f"Rounding this to the nearest 10 results in: {rounded_ferrite_level}")

if __name__ == "__main__":
    calculate_ferrite_level()