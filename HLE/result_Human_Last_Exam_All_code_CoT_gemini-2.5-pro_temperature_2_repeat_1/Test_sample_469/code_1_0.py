import math

def solve_color_problem():
    """
    This function determines the color of a Bromophenol blue solution
    at a specific pH when viewed through different path lengths.
    """
    # Step 1: Define the known values.
    # The pKa of Bromophenol Blue is ~4.1. Its acidic form is yellow, and its basic form is blue.
    ph = 4.21
    pka = 4.1
    path_length_thin_mm = 1
    path_length_thick_cm = 10
    
    # Convert all path lengths to cm for comparison.
    path_length_thin_cm = path_length_thin_mm / 10.0
    
    print("--- Analysis of Bromophenol Blue Solution Color ---")
    
    # Step 2: Determine the intrinsic color using the Henderson-Hasselbalch equation.
    print("\n1. Determining the solution's inherent color using pH and pKa.")
    print("The Henderson-Hasselbalch equation is: pH = pKa + log([Base Form]/[Acid Form])")
    print("For Bromophenol Blue, the Acid Form is yellow and the Base Form is blue.")
    print("\nPlugging in the numbers for the equation to find the ratio of forms:")
    print(f"{ph} = {pka} + log([Blue Form]/[Yellow Form])")
    
    log_ratio = ph - pka
    print(f"log([Blue Form]/[Yellow Form]) = {ph} - {pka} = {log_ratio:.2f}")
    
    ratio = math.pow(10, log_ratio)
    print(f"Ratio ([Blue Form]/[Yellow Form]) = 10^{log_ratio:.2f} = {ratio:.2f}")
    
    print("\nSince the ratio is close to 1, both the blue and yellow forms are present in significant amounts.")
    print("A mixture of yellow and blue light results in a perceived green color.")

    # Step 3: Apply the Beer-Lambert Law to account for path length.
    print("\n2. Considering the effect of path length (Beer-Lambert Law).")
    print("The intensity of perceived color is proportional to the path length of light through the solution.")
    
    print(f"\n- Thin side path length: {path_length_thin_cm} cm")
    print("This short path length means light absorption will be low, resulting in a pale or light color.")
    print("Resulting color: Light Green")
    
    print(f"\n- Thick side path length: {path_length_thick_cm} cm")
    print(f"This path is {path_length_thick_cm / path_length_thin_cm:.0f} times longer than the thin side.")
    print("This long path length means light absorption will be high, resulting in a deep, saturated color.")
    print("Resulting color: Green")
    
    print("\n--- Final Conclusion ---")
    print("Thin: light green, Thick: green")

solve_color_problem()
<<<C>>>