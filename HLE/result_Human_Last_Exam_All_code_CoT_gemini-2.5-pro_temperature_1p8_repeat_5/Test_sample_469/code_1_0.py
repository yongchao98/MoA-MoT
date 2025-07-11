import math

def analyze_solution_color():
    """
    Analyzes the color of a Bromophenol blue solution based on pH and viewing path length.
    """
    # Step 1: Define the problem parameters
    pKa_bromophenol_blue = 4.1
    pH_solution = 4.21
    path_length_thin_mm = 1
    path_length_thick_cm = 10

    # Convert all path lengths to cm for a consistent comparison
    path_length_thin_cm = path_length_thin_mm / 10.0

    # Step 2: Explain the intrinsic color using the Henderson-Hasselbalch equation
    # pH = pKa + log([Blue]/[Yellow])
    ratio_blue_to_yellow = 10**(pH_solution - pKa_bromophenol_blue)

    print("--- Analysis of Intrinsic Color ---")
    print(f"The pH of the solution is {pH_solution}.")
    print(f"The pKa of Bromophenol Blue is approximately {pKa_bromophenol_blue}.")
    print(f"The ratio of the blue form to the yellow form is 10^({pH_solution} - {pKa_bromophenol_blue}) = {ratio_blue_to_yellow:.2f}.")
    print("Since there are significant amounts of both yellow and blue forms, the solution appears green.\n")

    # Step 3: Explain the effect of path length using the Beer-Lambert Law
    absorbance_ratio = path_length_thick_cm / path_length_thin_cm

    print("--- Analysis of Color Intensity vs. Path Length ---")
    print(f"Path length of the thin side = {path_length_thin_cm} cm.")
    print(f"Path length of the thick side = {path_length_thick_cm} cm.")
    print("According to the Beer-Lambert Law, absorbance is proportional to path length.")
    print(f"The absorbance through the thick side is {path_length_thick_cm} cm / {path_length_thin_cm} cm = {absorbance_ratio:.0f} times greater than through the thin side.\n")

    # Step 4: Final Conclusion
    print("--- Conclusion ---")
    print("A shorter path length (thin side) results in lower absorbance, making the solution appear light green.")
    print("A longer path length (thick side) results in much higher absorbance, making the solution appear as a more saturated green.")
    print("\nFinal Answer: Thin: light green, Thick: green.")

analyze_solution_color()
<<<C>>>