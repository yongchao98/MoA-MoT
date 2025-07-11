import math

def calculate_brightness_drop():
    """
    Calculates the brightness drop of the brown dwarf during the transit of Planet 2.
    
    The calculation follows these steps:
    1. Determine the ratio of orbital radii (beta = r1/r2) based on the properties
       of the parabolic and circular orbits, leading to beta = 1/2.
    2. Use the geometric relations from the problem statement to find the ratio
       of Planet 2's radius to the brown dwarf's radius (R2/RB).
       R2/RB = (1 - beta) / (4 * beta)
    3. Calculate the squared ratio of radii, which corresponds to the ratio of areas.
    4. Use the magnitude formula to find the brightness drop in bolometric magnitudes.
       delta_m = -2.5 * log10(1 - (R2/RB)^2)
    """
    
    # Step 1: Ratio of orbital radii from orbital dynamics
    # This value is derived in the text explanation.
    beta = 0.5
    
    # Step 2: Ratio of Planet 2's radius to the Brown Dwarf's radius
    ratio_R2_RB = (1 - beta) / (4 * beta)
    
    # Step 3: Ratio of the projected areas
    ratio_areas = ratio_R2_RB ** 2
    
    # Step 4: Flux ratio during transit
    flux_ratio = 1 - ratio_areas
    
    # Step 5: Brightness drop in magnitudes
    delta_m = -2.5 * math.log10(flux_ratio)

    print("The ratio of the radius of Planet 2 to the radius of the brown dwarf (R2/RB) is: {:.2f}".format(ratio_R2_RB))
    print("The ratio of the projected areas (R2/RB)^2 is: {:.4f}".format(ratio_areas))
    print("\nThe brightness drop is calculated using the formula: Δm = -2.5 * log10(1 - (R2/RB)^2)")
    print("Plugging in the numbers:")
    print(f"Δm = -2.5 * log10(1 - {ratio_areas:.4f})")
    print(f"Δm = -2.5 * log10({flux_ratio:.4f})")

    # Print the final answer
    print(f"\nThe calculated brightness drop is: {delta_m:.5f} mag")
    print(f"The brightness drop rounded to three decimal places is: {delta_m:.3f} mag")
    
    final_answer = round(delta_m, 3)
    # The final answer needs to be enclosed in <<<>>>
    print(f"\n<<<{final_answer}>>>")

calculate_brightness_drop()