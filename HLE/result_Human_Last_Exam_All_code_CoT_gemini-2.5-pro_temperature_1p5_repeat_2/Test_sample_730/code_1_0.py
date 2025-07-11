import math

def solve():
    """
    This function lays out the calculation steps for the parrot.
    The goal is to estimate the mass of a spherical rock using fractions
    with integers no larger than 10.
    """

    # Explain the plan and the values being used.
    print("Yes, the mass can be estimated with an error of less than 10%.")
    print("Here is the calculation with the smallest possible integers.")
    print("-" * 20)
    print("The formula is: Mass = Density * (4/3) * pi * Radius^3")
    print("Given values: Density = 9/10, Radius = 1/2.")
    print("\nTo keep the numbers small, we simplify parts of the formula first:")
    
    # Step 1: Simplify Density * (4/3)
    print("Step 1: Combine density and the formula constant.")
    print("(9 / 10) * (4 / 3) = 36 / 30 = 6 / 5")

    # Step 2: Calculate Radius^3
    print("\nStep 2: Calculate the cube of the radius.")
    print("(1 / 2)^3 = (1 / 2) * (1 / 2) * (1 / 2) = 1 / 8")

    # Step 3: Choose an approximation for pi
    print("\nStep 3: Approximate pi with a simple fraction. We'll use 3/1.")
    
    # Final Calculation
    print("\nStep 4: Perform the final multiplication with the simplified parts.")
    print("Mass = (6 / 5) * (3 / 1) * (1 / 8)")
    print("     = (6 * 3 * 1) / (5 * 1 * 8)")
    print("     = 18 / 40")
    print("     = 9 / 20 kg")

    # Verify the largest integer used in the final step's expression
    # The calculation is (6/5) * (3/1) * (1/8). The integers are 6, 5, 3, 1, 8.
    # The largest integer is 8.

solve()