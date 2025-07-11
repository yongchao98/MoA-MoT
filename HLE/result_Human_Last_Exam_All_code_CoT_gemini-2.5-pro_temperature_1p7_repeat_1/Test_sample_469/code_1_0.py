import math

def solve_color_problem():
    """
    Determines the color of a Bromophenol blue solution in a rectangular prism
    viewed through its thin and thick sides.
    """
    # Step 1: Define constants and inputs
    pKa_bromophenol_blue = 4.1
    pH_solution = 4.21
    concentration_mg_ml = 13
    
    # Prism dimensions
    side1_mm = 1
    side2_cm = 10
    side3_cm = 10

    # Step 2: Determine the intrinsic color of the solution
    # The color of a pH indicator depends on the pH relative to its pKa.
    # Bromophenol blue is yellow below its pKa and blue above it.
    # When the pH is very close to the pKa, a mixture of both forms exists.
    # The mix of yellow and blue light appears green to the human eye.
    
    # Check if pH is in the transition range (approximated as pKa +/- 1)
    if pKa_bromophenol_blue - 1 < pH_solution < pKa_bromophenol_blue + 1:
        # Since pH 4.21 is very close to pKa 4.1, the color is green.
        intrinsic_color = "green"
    elif pH_solution <= pKa_bromophenol_blue - 1:
        intrinsic_color = "yellow"
    else: # pH_solution >= pKa_bromophenol_blue + 1
        intrinsic_color = "blue"

    # Step 3: Identify the path lengths
    # Convert all dimensions to a consistent unit (cm).
    side1_cm = side1_mm / 10.0
    
    # The path lengths for light are the side lengths of the prism.
    thin_path_cm = min(side1_cm, side2_cm, side3_cm)
    thick_path_cm = max(side1_cm, side2_cm, side3_cm)

    # Step 4: Apply the Beer-Lambert Law (A = εbc)
    # A = Absorbance, ε = molar absorptivity, b = path length, c = concentration
    # Absorbance is directly proportional to the path length (b).
    # A greater path length leads to higher absorbance and a more intense color.
    
    # We will show the relationship using the ratio of the path lengths.
    # Absorbance_thick / Absorbance_thin = (ε * b_thick * c) / (ε * b_thin * c) = b_thick / b_thin
    absorbance_ratio = thick_path_cm / thin_path_cm

    print(f"The intrinsic color of the solution at pH {pH_solution} is {intrinsic_color}.")
    print("-" * 30)
    print("The Beer-Lambert Law states that Absorbance is proportional to path length.")
    print(f"Path length of the 'thin' side = {thin_path_cm} cm")
    print(f"Path length of the 'thick' side = {thick_path_cm} cm")
    print("\nThe final equation for the ratio of absorbances is:")
    print(f"A_thick / A_thin = b_thick / b_thin = {thick_path_cm} cm / {thin_path_cm} cm = {absorbance_ratio:.0f}")
    print("-" * 30)
    
    # Step 5: Determine final color appearances
    # The shorter path length will result in a lighter, less saturated color.
    # The longer path length will result in a deeper, more saturated color.
    color_thin_side = f"light {intrinsic_color}"
    color_thick_side = f"{intrinsic_color}"
    
    print("Conclusion:")
    print(f"When viewed through the thin side, the color will be: {color_thin_side}")
    print(f"When viewed through the thick side, the color will be: {color_thick_side}")


solve_color_problem()