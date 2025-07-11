import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    based on initial conditions and forces.
    """
    # Define the physical constants and parameters of the problem
    # Gravitational Constant (m^3 kg^-1 s^-2)
    G = 6.67430e-11
    # Mass of the asteroid A (kg)
    m = 5.972e20
    # Mass of the spaceship B (kg)
    M = 10000.0
    # Initial distance between A and B (m)
    l_0 = 1.0e6
    # Initial speed of the spaceship (m/s)
    v_0 = 100.0
    # Additional applied force on the spaceship (N)
    F = 100.0

    print("Problem Parameters:")
    print(f"  Gravitational Constant (G): {G:.4e} m^3 kg^-1 s^-2")
    print(f"  Asteroid Mass (m): {m:.4e} kg")
    print(f"  Spaceship Mass (M): {M:.4e} kg")
    print(f"  Initial Distance (l_0): {l_0:.4e} m")
    print(f"  Initial Speed (v_0): {v_0:.4e} m/s")
    print(f"  Applied Force (F): {F:.4e} N\n")
    
    # Based on the work-energy theorem, we derive a quadratic equation for l_max:
    # a * (l_max)^2 + b * l_max + c = 0
    
    # Calculate the coefficients a, b, and c
    a = F
    b = -(F * l_0 + (G * M * m) / l_0 - 0.5 * M * v_0**2)
    c = G * M * m

    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("The conditions (e.g., high initial speed or strong applied force) are such that")
        print("the spaceship has enough energy to escape to infinity. No maximum distance is reached.")
        return

    # Use the quadratic formula to find the two possible distances (roots)
    sqrt_discriminant = math.sqrt(discriminant)
    root1 = (-b + sqrt_discriminant) / (2 * a)
    root2 = (-b - sqrt_discriminant) / (2 * a)
    
    # The maximum distance is the first root encountered that is greater than the initial distance
    # Since the spaceship is moving away from l_0, we look for the smallest root > l_0.
    l_max = -1
    if root1 >= l_0 and root2 >= l_0:
        l_max = min(root1, root2)
    elif root1 >= l_0:
        l_max = root1
    elif root2 >= l_0:
        l_max = root2

    print("The derived quadratic equation for the maximum distance (l_max) is:")
    print(f"({a:.4e}) * l_max^2 + ({b:.4e}) * l_max + ({c:.4e}) = 0\n")

    if l_max != -1:
        print(f"The two solutions for the distance are {root1:.2f} m and {root2:.2f} m.")
        print(f"Starting from l_0 = {l_0:.2f} m, the first turning point reached is the maximum distance.")
        print(f"\nThe calculated maximum distance is {l_max:.2f} meters.")
        # Final answer format for the platform
        print(f"\n<<<{l_max:.1f}>>>")
    else:
        print("Could not find a physical solution. The spaceship might not be able to move away from the asteroid under these conditions.")

if __name__ == "__main__":
    calculate_max_distance()