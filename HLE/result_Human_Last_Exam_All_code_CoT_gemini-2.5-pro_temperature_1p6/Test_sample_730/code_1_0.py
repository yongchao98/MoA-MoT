import math

def parrot_calculation():
    """
    This function outlines the calculation for the parrot to estimate the rock's mass.
    """
    print("Yes, the parrot can estimate the mass. Here is the recommended calculation:")
    print("-" * 50)

    # 1. State the formula and given values
    print("The formula is: Mass = (4/3) * π * (radius)^3 * density")
    print("\nWe are given:")
    print("radius (r) = 1/2 cm")
    print("density (ρ) = 9/10 kg/cm^3")

    # 2. State the approximation for Pi
    print("\nTo keep all integers small and the error below 10%, we will approximate π ≈ 10/3.")
    # Actual mass = (4/3) * pi * (1/2)**3 * (9/10) = (3/20)*pi approx 0.4712 kg
    # Estimated mass = (4/3) * (10/3) * (1/2)**3 * (9/10) = 1/2 = 0.5 kg
    # Error = abs(0.5 - 0.4712) / 0.4712 is approx 6.1%, which is less than 10%.

    # 3. Show the full calculation with all numbers
    print("\nThe full calculation with all numbers substituted is:")
    print("Mass ≈ (4/3) * (10/3) * (1/2)^3 * (9/10)")

    print("\nThe integers involved in this calculation are: 4, 3, 10, 3, 1, 2, 9, 10")

    # 4. Show the step-by-step simplification
    print("\nHere is a simple way to calculate this by cancelling terms:")
    print("Step 1: Write out all the fractions multiplied together.")
    print("  Mass ≈ (4/3) * (10/3) * (1/8) * (9/10)")
    print("\nStep 2: Group the numerators and denominators.")
    print("  Mass ≈ (4 * 10 * 1 * 9) / (3 * 3 * 8 * 10)")
    print("\nStep 3: Simplify by cancelling out common terms.")
    print("  The (3 * 3) in the bottom cancels with the 9 on top.")
    print("  The 10 on the bottom cancels with the 10 on top.")
    print("  The 4 on top and 8 on the bottom simplify to 1/2.")
    print("  Mass ≈ 4 / 8")
    print("\nStep 4: Final Result.")
    print("  Mass ≈ 1/2 kg")
    print("-" * 50)

parrot_calculation()
