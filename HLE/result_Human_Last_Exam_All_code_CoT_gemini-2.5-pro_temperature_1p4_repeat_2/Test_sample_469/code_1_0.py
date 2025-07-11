import math

def solve_color_problem():
    """
    Determines the color of a Bromophenol blue solution at a specific pH
    when viewed through different path lengths.
    """
    # --- Step 1: Define the known values ---
    pKa_bromophenol_blue = 4.0
    pH_solution = 4.21
    path_length_thin_mm = 1.0
    path_length_thick_cm = 10.0

    # --- Step 2: Determine the inherent color of the solution ---
    # We use the Henderson-Hasselbalch equation to find the ratio of the indicator's forms.
    # For Bromophenol blue, the Acid form (HIn) is YELLOW and the Base form (In-) is BLUE.
    # Equation: pH = pKa + log([Base]/[Acid])
    
    # Rearranging the equation: log([Base]/[Acid]) = pH - pKa
    log_ratio = pH_solution - pKa_bromophenol_blue

    # Ratio = 10^(log_ratio)
    ratio_base_to_acid = 10**log_ratio

    print("--- Part 1: Determining the Solution's Inherent Color ---")
    print("The Henderson-Hasselbalch equation is: pH = pKa + log([Base Form]/[Acid Form])")
    print(f"For this solution, the equation with the given values is:")
    print(f"log([Blue Form]/[Yellow Form]) = {pH_solution} - {pKa_bromophenol_blue}")
    log_ratio_val = pH_solution - pKa_bromophenol_blue
    print(f"log([Blue Form]/[Yellow Form]) = {log_ratio_val:.2f}")
    ratio_val = 10**log_ratio_val
    print(f"The ratio of [Blue Form]/[Yellow Form] is 10^{log_ratio_val:.2f} which is approximately {ratio_val:.2f}.")
    print("\nBecause there are significant amounts of both the YELLOW acid form and the BLUE base form,")
    print("the transmitted light will appear GREEN to the human eye.")

    # --- Step 3: Analyze the effect of path length on color intensity ---
    # Beer-Lambert Law: Absorbance (A) is proportional to path length (b).
    
    # Convert path lengths to a consistent unit (cm)
    path_length_thin_cm = path_length_thin_mm / 10.0
    
    print("\n--- Part 2: Analyzing the Effect of Path Length ---")
    print("According to the Beer-Lambert Law, the intensity of the perceived color depends on the path length.")
    print(f"Thin side path length = {path_length_thin_mm} mm = {path_length_thin_cm} cm")
    print(f"Thick side path length = {path_length_thick_cm} cm")
    
    absorbance_ratio = path_length_thick_cm / path_length_thin_cm
    print(f"The absorbance through the thick side is {absorbance_ratio:.0f} times greater than through the thin side.")
    
    # --- Step 4: Conclusion ---
    print("\n--- Final Conclusion ---")
    print("Viewing through the THIN side (short path, low absorbance) results in a LIGHT GREEN color.")
    print("Viewing through the THICK side (long path, high absorbance) results in a more intense GREEN color.")

solve_color_problem()
<<<C>>>