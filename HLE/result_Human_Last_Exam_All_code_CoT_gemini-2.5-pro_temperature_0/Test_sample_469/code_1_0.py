import math

def solve_color_problem():
    """
    Determines the color of a Bromophenol blue solution based on pH and path length.
    """
    # --- Step 1: Define constants and properties ---
    # Properties of Bromophenol Blue
    pKa = 4.1
    color_acidic = "yellow"  # Color of the HIn form (at low pH)
    color_basic = "blue"     # Color of the In- form (at high pH)
    transition_color = "green" # The color when yellow and blue forms are mixed

    # Solution and prism properties
    pH = 4.21
    path_length_thin_mm = 1
    path_length_thick_cm = 10

    # --- Step 2: Convert units for consistency ---
    path_length_thin_cm = path_length_thin_mm / 10.0

    # --- Step 3: Determine the ratio of indicator forms using the Henderson-Hasselbalch equation ---
    # The equation is: pH = pKa + log10([Base Form]/[Acid Form])
    # Rearranging for the ratio: [Base Form]/[Acid Form] = 10^(pH - pKa)
    # [Base Form] is the concentration of the basic (blue) form.
    # [Acid Form] is the concentration of the acidic (yellow) form.

    exponent = pH - pKa
    ratio_base_to_acid = 10**exponent

    # --- Step 4: Print the analysis and conclusion ---
    print("Analysis of Bromophenol Blue Solution Color")
    print("="*50)

    print(f"1. Determine the intrinsic color of the solution at pH {pH}:")
    print(f"   - The pKa of Bromophenol Blue is {pKa}.")
    print(f"   - The acidic form is {color_acidic} and the basic form is {color_basic}.")
    print("\n   Using the Henderson-Hasselbalch equation to find the ratio of the forms:")
    print(f"   Ratio ([Blue Form] / [Yellow Form]) = 10^(pH - pKa)")
    # Outputting each number in the final equation
    print(f"   Ratio = 10^({pH} - {pKa})")
    print(f"   Ratio = 10^({exponent:.2f}) ≈ {ratio_base_to_acid:.2f}")
    print("\n   Since the pH is in the indicator's transition range, both the yellow")
    print(f"   and blue forms are present. The mixture of these colors appears {transition_color}.")
    print("-" * 50)

    print("2. Consider the effect of viewing path length (Beer-Lambert Law):")
    print(f"   - Thin side path length: {path_length_thin_mm} mm = {path_length_thin_cm} cm")
    print(f"   - Thick side path length: {path_length_thick_cm} cm")
    print("\n   According to the Beer-Lambert Law, absorbance is proportional to path length.")
    print(f"   A longer path length results in higher absorbance and a more intense color.")
    print("-" * 50)

    print("3. Conclusion:")
    print(f"   - Through the short path ({path_length_thin_cm} cm), the color will be a light {transition_color}.")
    print(f"   - Through the long path ({path_length_thick_cm} cm), the color will be a more intense, deep {transition_color}.")
    print("\n   This corresponds to the answer: Thin: light green, Thick: green")

solve_color_problem()