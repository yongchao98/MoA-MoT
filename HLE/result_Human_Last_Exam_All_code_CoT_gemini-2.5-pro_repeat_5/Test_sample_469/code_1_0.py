import math

def solve_color_problem():
    """
    This function determines the color of a Bromophenol blue solution
    at a specific pH when viewed through different path lengths.
    """
    # Given values
    ph = 4.21
    pka_bromophenol_blue = 4.0
    path_length_thin_mm = 1
    path_length_thick_cm = 10

    # Convert path lengths to a consistent unit (cm)
    path_length_thin_cm = path_length_thin_mm / 10.0

    print("Step 1: Determine the intrinsic color of the solution at pH {}.".format(ph))
    print("The pKa of Bromophenol blue is approximately {}.".format(pka_bromophenol_blue))
    print("Its transition range is from yellow (acidic) to blue (basic).")
    print("Since pH {:.2f} is close to the pKa, the solution will be a mix of both forms.".format(ph))
    print("A mixture of yellow and blue appears green.")

    print("\nWe can use the Henderson-Hasselbalch equation to confirm this.")
    print("The equation is: pH = pKa + log([Base Form] / [Acid Form])")
    
    # Calculate the ratio of the base form (blue) to the acid form (yellow)
    ratio_base_to_acid = math.pow(10, ph - pka_bromophenol_blue)
    
    print("\nFinal Equation with values:")
    print("{:.2f} = {:.1f} + log([Blue Form] / [Yellow Form])".format(ph, pka_bromophenol_blue))
    print("The ratio of the blue form to the yellow form is {:.2f}:1.".format(ratio_base_to_acid))
    print("This significant mixture confirms the color is green.")

    print("\nStep 2: Analyze the effect of path length (Beer-Lambert Law).")
    print("Absorbance is proportional to path length. A longer path means a more intense color.")
    print("- Thin side path length: {} cm".format(path_length_thin_cm))
    print("- Thick side path length: {} cm".format(path_length_thick_cm))
    print("The thick side has a path length {} times longer than the thin side.".format(path_length_thick_cm / path_length_thin_cm))

    print("\nConclusion:")
    print("- Through the thin side (low absorbance), the color will be light green.")
    print("- Through the thick side (high absorbance), the color will be a more intense green.")

solve_color_problem()
<<<C>>>