import math
from fractions import Fraction

def calculate_stress_ratio():
    """
    Calculates the ratio of maximum tangential stress to internal pressure
    for a thick-walled cylinder where the outer radius is twice the inner radius.
    """
    print("In elasticity theory, the maximum tangential (hoop) stress in a thick-walled cylinder under internal pressure (Pi) occurs at the inner radius.")
    print("The formula for this maximum stress (σ_t_max) is derived from Lamé's equations:")
    print("σ_t_max = Pi * (r_o^2 + r_i^2) / (r_o^2 - r_i^2)")
    print("\nWe need to find the ratio of σ_t_max to Pi:")
    print("Ratio = σ_t_max / Pi = (r_o^2 + r_i^2) / (r_o^2 - r_i^2)")

    # Define the relationship between the radii.
    # Let's set the inner radius (r_i) to 1 for simplicity.
    r_i = 1.0
    # The outer radius (r_o) is twice the inner radius.
    r_o = 2.0 * r_i

    print(f"\nGiven that the outer radius (r_o) is twice the inner radius (r_i), we have r_o = {int(r_o)} and r_i = {int(r_i)}.")

    # Calculate the terms in the ratio formula
    r_o_sq = r_o**2
    r_i_sq = r_i**2

    # Calculate numerator and denominator
    numerator = r_o_sq + r_i_sq
    denominator = r_o_sq - r_i_sq

    # Calculate the final ratio
    ratio = numerator / denominator

    # Display the final equation with the calculated numbers
    print("\nSubstituting these values into the ratio formula:")
    print(f"Ratio = ({r_o_sq} + {r_i_sq}) / ({r_o_sq} - {r_i_sq})")
    print(f"Ratio = {numerator} / {denominator}")

    # Display the result as a fraction and a decimal
    ratio_fraction = Fraction(numerator, denominator)
    print(f"\nThe calculated ratio as a fraction is: {ratio_fraction.numerator}/{ratio_fraction.denominator}")
    print(f"The calculated ratio as a decimal is: {ratio:.3f}")

calculate_stress_ratio()