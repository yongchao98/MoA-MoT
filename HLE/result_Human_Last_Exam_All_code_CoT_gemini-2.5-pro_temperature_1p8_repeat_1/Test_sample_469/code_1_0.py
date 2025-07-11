import math

def solve_color_problem():
    """
    This script determines the color of a Bromophenol blue solution at a specific pH
    when viewed through different path lengths.
    """
    # --- Given Parameters ---
    # The pH of the aqueous Bromophenol blue solution
    ph_solution = 4.21
    # The approximate pKa of Bromophenol blue (the pH at which it is in transition)
    pka_bromophenol_blue = 4.1
    # Path length of the thin side of the prism
    path_thin_mm = 1
    # Path length of the thick side of the prism
    path_thick_cm = 10.0

    # --- Unit Conversion ---
    # Convert millimeters to centimeters to have consistent units
    path_thin_cm = path_thin_mm / 10.0

    # --- Analysis ---
    print("Step 1: Determine the intrinsic color of the solution using the Henderson-Hasselbalch equation.")
    print("The Henderson-Hasselbalch equation is: pH = pKa + log([Blue Form] / [Yellow Form])")
    print("We can find the ratio of the blue form to the yellow form.")
    
    # Calculate the ratio [Blue Form]/[Yellow Form] = 10^(pH - pKa)
    exponent = ph_solution - pka_bromophenol_blue
    ratio_blue_to_yellow = math.pow(10, exponent)

    print(f"The final equation for the ratio is: [Blue Form] / [Yellow Form] = 10^({ph_solution} - {pka_bromophenol_blue})")
    print(f"The ratio of the blue form to the yellow form is ~{ratio_blue_to_yellow:.2f}:1.")
    print("Since the ratio is close to 1, both yellow and blue forms are present.")
    print("The mix of yellow and blue results in a GREEN color.\n")

    print("Step 2: Analyze the effect of path length on color intensity using the Beer-Lambert Law (A = εbc).")
    print("The law states that Absorbance (color intensity) is proportional to the path length.")
    
    print(f"Path length of the thin side: {path_thin_mm} mm ({path_thin_cm} cm)")
    print(f"Path length of the thick side: {path_thick_cm} cm")
    
    # Calculate the ratio of the path lengths
    path_ratio = path_thick_cm / path_thin_cm
    print(f"The thick side's path is {int(path_ratio)} times longer than the thin side's.")
    print("A shorter path length (thin side) means lower absorbance and a lighter color.")
    print("A longer path length (thick side) means higher absorbance and a more intense, saturated color.\n")

    print("--- Conclusion ---")
    print("Thin side view: A GREEN color, but faint due to the short path length -> Light green.")
    print("Thick side view: A GREEN color, which is intense due to the long path length -> Green.")

solve_color_problem()
<<<C>>>