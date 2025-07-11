import math

def solve_color_problem():
    """
    This script determines the color of a Bromophenol blue solution at a specific pH
    when viewed through different path lengths of a rectangular prism.
    """
    # --- Step 1: Define physical and chemical properties ---
    ph_solution = 4.21
    # Bromophenol blue has a pKa of ~4.0 and transitions from yellow (acid) to blue (base)
    # in the pH range of 3.0 to 4.6. The mix of yellow and blue appears green.
    pka = 4.0
    color_acid = "yellow"
    color_base = "blue"
    color_transition = "green"
    
    # Prism dimensions and path lengths
    side1_mm = 1
    side2_cm = 10
    side3_cm = 10
    
    # Convert all dimensions to cm for consistency
    path_thin_cm = side1_mm / 10.0
    path_thick_cm = side2_cm

    # --- Step 2: Determine the intrinsic color based on pH ---
    print("Step 1: Determining the intrinsic color of the solution.")
    print(f"The pH of the solution is {ph_solution}.")
    print(f"The pKa of Bromophenol blue is approximately {pka}.")
    
    # Using the Henderson-Hasselbalch equation: pH = pKa + log10([base]/[acid])
    # The ratio of the basic (blue) form to the acidic (yellow) form can be calculated.
    ratio = 10**(ph_solution - pka)
    
    print(f"The ratio of the basic '{color_base}' form to the acidic '{color_acid}' form is 10^({ph_solution} - {pka:.1f}) = {ratio:.2f}.")
    print(f"Since the pH of {ph_solution} is within the indicator's transition range (3.0-4.6), both the yellow and blue forms are present.")
    print(f"The mixture of {color_acid} and {color_base} light results in a {color_transition} color.\n")

    # --- Step 3: Analyze the effect of path length (Beer-Lambert Law) ---
    print("Step 2: Analyzing the effect of the viewing path length.")
    print(f"The 'thin side' corresponds to a path length of {path_thin_cm} cm.")
    print(f"The 'thick side' corresponds to a path length of {path_thick_cm} cm.")
    print("The Beer-Lambert Law states that light absorbance is proportional to the path length.")
    print("A shorter path length results in lower absorbance, making the color appear lighter.")
    print("A longer path length results in higher absorbance, making the color appear more intense and saturated.\n")

    # --- Step 4: Conclude the final answer ---
    print("Step 3: Determining the final perceived colors.")
    color_thin_side = f"light {color_transition}"
    color_thick_side = f"{color_transition}"
    
    print(f"Through the thin side ({path_thin_cm} cm), the color will be: {color_thin_side}.")
    print(f"Through the thick side ({path_thick_cm} cm), the color will be: {color_thick_side}.")

solve_color_problem()
<<<C>>>