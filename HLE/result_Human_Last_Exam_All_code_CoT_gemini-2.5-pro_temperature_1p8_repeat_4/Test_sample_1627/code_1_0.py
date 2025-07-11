import fractions

def solve_cylinder_stress_ratio():
    """
    Calculates the ratio of maximum tangential stress to internal pressure
    for a thick-walled cylinder where the outer radius is twice the inner radius.
    """
    
    # The formula for the ratio is: (r_o^2 + r_i^2) / (r_o^2 - r_i^2)
    # where r_o is the outer radius and r_i is the inner radius.

    # Given condition: r_o = 2 * r_i.
    # We can set r_i = 1 for calculation, which means r_o = 2.
    # The actual radius value cancels out, so this simplification is valid.
    r_i = 1
    r_o = 2
    
    # Calculate the numerator and denominator of the ratio formula.
    numerator = r_o**2 + r_i**2
    denominator = r_o**2 - r_i**2
    
    # Display the problem setup and formula
    print("In elasticity theory, for a thick-walled cylinder under internal pressure,")
    print("the ratio of maximum tangential stress to internal pressure is given by:")
    print("Ratio = (r_o^2 + r_i^2) / (r_o^2 - r_i^2)")
    print("\nGiven r_o = 2 * r_i, we can set r_i=1 and r_o=2 for our calculation.")
    print("\nSubstituting these values into the formula:")

    # Print the equation with numbers plugged in
    print(f"Ratio = ({r_o}^2 + {r_i}^2) / ({r_o}^2 - {r_i}^2)")
    
    # Print the intermediate step with squared values
    print(f"Ratio = ({r_o**2} + {r_i**2}) / ({r_o**2} - {r_i**2})")

    # Print the final fraction before division
    print(f"Ratio = {numerator} / {denominator}")
    
    # Represent the final answer as a fraction
    final_ratio_fraction = fractions.Fraction(numerator, denominator)
    
    print(f"\nTherefore, the final ratio is {final_ratio_fraction}.")
    print(f"As a decimal, this is approximately {float(final_ratio_fraction):.3f}.")

solve_cylinder_stress_ratio()
<<<E>>>