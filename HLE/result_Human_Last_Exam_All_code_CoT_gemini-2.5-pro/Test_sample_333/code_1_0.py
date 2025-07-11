import sys

def solve_projectile_explosion():
    """
    Calculates the maximum landing distance of a fragment from an exploding projectile.
    """
    # Given horizontal distance to the highest point of elevation
    I = 500  # in meters

    # --- Explanation and Step-by-step Calculation ---

    print("This problem can be solved using the concept of the center of mass (CM).")
    print("The explosion is an internal force, so the trajectory of the CM is unaffected.")
    print("-" * 50)

    # Step 1: Determine the landing point of the Center of Mass
    # The full range of the original projectile would be 2 * I.
    # This is where the CM of the fragments will land.
    x_cm = 2 * I
    print(f"1. The projectile's highest point is at I = {I} m.")
    print(f"2. The Center of Mass (CM) continues on the original path and lands at x_cm = 2 * I.")
    print(f"   x_cm = 2 * {I} = {x_cm} m.")
    print("-" * 50)

    # Step 2: Use the Center of Mass formula for two equal masses
    # For two equal masses, x_cm = (x1 + x2) / 2
    # We want to find x2. Let's rearrange the formula:
    # 2 * x_cm = x1 + x2
    # x2 = 2 * x_cm - x1
    print("3. The CM position for two equal fragments is: x_cm = (x1 + x2) / 2.")
    print("   Where x1 and x2 are the landing spots of the fragments.")
    print("   Rearranging for x2 gives: x2 = 2 * x_cm - x1.")
    print("-" * 50)

    # Step 3: Define the landing position of the first fragment
    # To get the MAXIMUM distance for the second fragment, the first fragment
    # must land as close to the gun as possible. We assume it lands at the origin.
    x1 = 0  # in meters
    print(f"4. One fragment falls 'near the gun'. For maximum range of the other fragment,")
    print(f"   we assume this fragment lands at the gun's location, so x1 = {x1} m.")
    print("-" * 50)

    # Step 4: Calculate the final position of the second fragment
    print("5. Now, we substitute the known values into the equation to find x2.")
    
    # Calculate final answer
    x2 = 2 * x_cm - x1

    # Print the final equation with all numbers
    print("\nFinal Equation:")
    print(f"x2 = 2 * ({x_cm}) - {x1}")
    print(f"x2 = {x2}")

    print(f"\nTherefore, the maximum distance from the gun you can safely be is {x2} meters.")
    
    # Final answer in the required format
    sys.stdout.write(f"\n<<<{x2}>>>")

solve_projectile_explosion()