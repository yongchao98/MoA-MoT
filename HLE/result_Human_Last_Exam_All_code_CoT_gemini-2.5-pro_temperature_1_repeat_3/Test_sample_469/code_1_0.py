import math

def solve_color_problem():
    """
    This script solves the Bromophenol blue color problem by following these steps:
    1. Define the chemical and physical constants from the problem.
    2. Use the Henderson-Hasselbalch equation to determine the ratio of the different colored forms of Bromophenol blue at the given pH.
    3. Determine the resulting color from the mixture of forms.
    4. Use the Beer-Lambert law concept to explain the difference in color intensity between the thin and thick sides of the prism.
    """

    # 1. Define constants
    pH = 4.21
    pKa_bromophenol_blue = 4.1
    # The acidic form (HIn) is yellow, the basic form (In-) is blue.

    path_length_thin_mm = 1
    path_length_thick_cm = 10

    # Convert path lengths to cm for consistency
    path_length_thin_cm = path_length_thin_mm / 10.0

    # 2. Use the Henderson-Hasselbalch equation to find the ratio of the colored forms.
    # The equation is: pH = pKa + log([In-]/[HIn]), where [In-] is the blue form and [HIn] is the yellow form.
    
    print("Step 1: Determine the ratio of the blue to yellow forms of Bromophenol Blue.")
    print("We use the Henderson-Hasselbalch equation: pH = pKa + log([Blue Form]/[Yellow Form])")
    print("Plugging in the given values:")
    print(f"{pH} = {pKa_bromophenol_blue} + log([Blue Form]/[Yellow Form])")
    
    log_ratio = pH - pKa_bromophenol_blue
    print(f"log([Blue Form]/[Yellow Form]) = {pH} - {pKa_bromophenol_blue} = {log_ratio:.2f}")

    ratio = math.pow(10, log_ratio)
    print(f"The ratio of [Blue Form]/[Yellow Form] is 10^{log_ratio:.2f} = {ratio:.2f}")
    print("-" * 20)

    # 3. Determine the resulting color from the mixture.
    print("Step 2: Determine the fundamental color of the solution.")
    print(f"The ratio of the blue form to the yellow form is approximately {ratio:.2f}.")
    print("Since the pH is very close to the pKa, there are significant amounts of both the yellow acidic form and the blue basic form in the solution.")
    print("When yellow and blue light are mixed, the resulting perceived color is green.")
    print("Therefore, the fundamental color of the solution is green.")
    print("-" * 20)

    # 4. Explain the effect of path length on color intensity.
    print("Step 3: Analyze the effect of the viewing path length.")
    print("The Beer-Lambert Law (A = εbc) states that absorbance (A) is directly proportional to the path length (b).")
    print(f"The thin side has a path length of {path_length_thin_cm} cm.")
    print(f"The thick side has a path length of {path_length_thick_cm} cm.")
    
    path_ratio = path_length_thick_cm / path_length_thin_cm
    print(f"The thick side is {path_ratio:.0f} times longer than the thin side.")
    
    print("\nConclusion:")
    print(" - Through the THIN side (short path length), the absorbance is lower, resulting in a less intense, or 'light green' color.")
    print(" - Through the THICK side (long path length), the absorbance is much higher, resulting in a more intense, darker 'green' color.")
    print("\nComparing this to the choices, the best description is 'Thin: light green, Thick: green'.")

solve_color_problem()
<<<C>>>