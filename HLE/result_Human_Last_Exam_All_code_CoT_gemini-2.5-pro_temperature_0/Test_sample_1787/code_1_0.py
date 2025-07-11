import math

def calculate_brightness_drop():
    """
    Calculates the brightness drop of a brown dwarf during a planetary transit based on the problem's parameters.
    """
    # Step 1 & 2: The relationship between the radius ratio and distance ratio is derived in the text.
    # R_2 / R_BD = (1 - x) / (4 * x), where x = d_1 / d_2.

    # Step 3: Determine the distance ratio x = d_1 / d_2.
    # Based on the orbital properties (P1 parabolic at pericenter, P2 circular),
    # we assume the specific angular momenta are equal.
    # h_1 = sqrt(2*G*M*d_1)
    # h_2 = sqrt(G*M*d_2)
    # h_1 = h_2  =>  2*d_1 = d_2  =>  x = d_1 / d_2 = 0.5
    x = 0.5
    print(f"The ratio of the orbital distances, x = d1/d2, is determined to be: {x}")

    # Step 4: Calculate the ratio of the radii of Planet 2 and the Brown Dwarf.
    ratio_R2_RBD = (1 - x) / (4 * x)
    print(f"The ratio of the radius of Planet 2 to the Brown Dwarf, R2/R_BD, is: {ratio_R2_RBD:.4f}")

    # Step 5: Calculate the brightness drop in magnitudes.
    # The formula is delta_m = -2.5 * log10(1 - (R2/R_BD)^2)
    flux_ratio = 1 - ratio_R2_RBD**2
    
    print("\nThe brightness drop is calculated using the formula: delta_m = -2.5 * log10(1 - (R2/R_BD)^2)")
    print(f"Plugging in the numbers: delta_m = -2.5 * log10(1 - ({ratio_R2_RBD:.4f})^2)")
    print(f"This simplifies to: delta_m = -2.5 * log10({flux_ratio:.4f})")

    delta_m = -2.5 * math.log10(flux_ratio)
    
    print(f"\nThe calculated brightness drop is {delta_m:.5f} magnitudes.")
    print(f"Rounded to three decimal places, the final answer is {delta_m:.3f} magnitudes.")
    
    # Final answer in the required format
    print(f"\n<<<The brightness drop is {delta_m:.3f} magnitudes.>>>")


calculate_brightness_drop()