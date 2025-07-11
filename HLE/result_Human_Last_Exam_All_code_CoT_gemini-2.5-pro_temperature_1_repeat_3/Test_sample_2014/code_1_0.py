import math

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two tandem aerofoils in ground effect.
    """
    print("This script calculates the lift ratio L1/L2 for two aerofoils in tandem near the ground.")
    
    # We can set the chord 'c' to 1.0, as the geometry is given in terms of c,
    # and it will cancel out in the final dimensionless ratio.
    c = 1.0
    print(f"Step 1: Define the geometry based on the chord c = {c}\n")

    # Given geometric parameters from the problem statement
    s = 0.5 * c
    h = 0.5 * c
    print(f"The separation between aerofoils, s = 1/2 * c = {s}")
    print(f"The ride height above the ground, h = 1/2 * c = {h}\n")

    # The lift ratio L1/L2 is determined by an interaction factor 'A'.
    # This factor 'A' is derived from vortex theory.
    print("Step 2: Calculate the interaction factor 'A'.")
    print("The equation for A is: A = (c/2) * [ 1/s - s / (s^2 + 4*h^2) ]\n")

    # Calculate the individual terms to show the process clearly
    term1 = 1 / s
    term2_num = s
    term2_den = s**2 + 4 * h**2
    term2 = term2_num / term2_den
    
    print("Calculating the parts of the equation for A:")
    print(f"1/s = 1/{s} = {term1}")
    print(f"s / (s^2 + 4*h^2) = {s} / ({s}^2 + 4*{h}^2) = {s} / ({s**2:.2f} + {4*h**2:.2f}) = {s} / {term2_den:.2f} = {term2:.3f}\n")

    # Calculate the final value of A
    A = (c / 2.0) * (term1 - term2)
    
    print("Substituting these values to find A:")
    print(f"A = ({c}/2) * [ {term1} - {term2:.3f} ]")
    print(f"A = {c/2} * [ {term1 - term2:.3f} ]")
    print(f"A = {A:.3f}\n")

    # The final step is to calculate the lift ratio using the value of A.
    print("Step 3: Calculate the lift ratio L1/L2 using A.")
    print("The equation for the lift ratio is: L1/L2 = (1 - A) / (1 + A)\n")

    # Calculate the numerator and denominator of the ratio
    ratio_num = 1.0 - A
    ratio_den = 1.0 + A
    lift_ratio = ratio_num / ratio_den

    print("Substituting the value of A:")
    print(f"L1/L2 = (1 - {A:.3f}) / (1 + {A:.3f})")
    print(f"L1/L2 = {ratio_num:.3f} / {ratio_den:.3f}")
    print(f"L1/L2 = {lift_ratio:.4f}\n")

    print("The final result is approximately 0.1111, which is equal to the fraction 1/9.")

calculate_lift_ratio()
<<<0.1111>>>