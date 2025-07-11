import math

def calculate_ferrite_level():
    """
    Calculates the approximate ferrite level for a stainless steel using the
    Kotecki and Siewert formula, an approximation of the WRC-1992 diagram.
    """
    # Input values
    cr_eq = 39.0
    ni_eq = 29.0

    # Step 1: Calculate the chromium contribution term (FN_Cr)
    # Formula: FN_Cr = 4.88 * Cr_eq - 49.9
    fn_cr = 4.88 * cr_eq - 49.9

    # Step 2: Calculate the nickel contribution term (FN_Ni)
    # Formula: FN_Ni = 0.05 * Ni_eq^2 + 0.52 * Ni_eq
    fn_ni = 0.05 * (ni_eq ** 2) + 0.52 * ni_eq

    # Step 3: Define constants for the final formula
    constant_term = 6.6
    divisor = 0.82
    
    # Step 4: Calculate the final Ferrite Level (FN)
    # Formula: FN = (FN_Cr - FN_Ni - 6.6) / 0.82
    ferrite_level = (fn_cr - fn_ni - constant_term) / divisor

    # Step 5: Round the result to the nearest 10
    # The rounding logic: round(value / 10) * 10
    rounded_ferrite_level = round(ferrite_level / 10) * 10

    # Print the process and the final equation with all numerical values
    print("To find the ferrite level, we use the Kotecki and Siewert formula with:")
    print(f"Chromium Equivalent (Cr_eq) = {cr_eq}")
    print(f"Nickel Equivalent (Ni_eq) = {ni_eq}\n")
    
    print("The final calculation is broken down as follows:")
    print(f"Ferrite Level = ((4.88 * {cr_eq} - 49.9) - (0.05 * {ni_eq}**2 + 0.52 * {ni_eq}) - {constant_term}) / {divisor}")
    print(f"Ferrite Level = ({fn_cr:.2f} - {fn_ni:.2f} - {constant_term}) / {divisor}")
    print(f"Ferrite Level = {(fn_cr - fn_ni - constant_term):.2f} / {divisor}")
    print(f"Ferrite Level = {ferrite_level:.2f}\n")
    
    print(f"The calculated ferrite level is approximately {ferrite_level:.2f}.")
    print(f"When rounded to the nearest 10, the result is: {int(rounded_ferrite_level)}")

if __name__ == '__main__':
    calculate_ferrite_level()