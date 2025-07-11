import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    given gravitational pull and a constant counteracting force.
    """
    # --- Problem Parameters ---
    # You can change these values to test different scenarios.

    # Mass of the asteroid A (in kg)
    m = 5.972e24  # e.g., mass of Earth
    # Mass of the spaceship B (in kg)
    M = 100000.0   # 100,000 kg
    # Initial distance between A and B (in meters)
    l_0 = 7.0e6      # e.g., 7,000 km from center of mass A
    # Initial speed of spaceship B away from A (in m/s)
    v_0 = 1000.0     # 1 km/s
    # Additional constant force on spaceship B (in Newtons)
    F = 5.0e5      # 500,000 N

    # Physical constant
    # Gravitational constant (m^3 kg^-1 s^-2)
    G = 6.67430e-11

    # --- Calculation ---

    # Based on the conservation of energy, we derive a quadratic equation
    # for the maximum distance l_max of the form:
    # a*(l_max)^2 + b*(l_max) + c = 0

    GmM = G * m * M
    
    # Coefficient 'a' of the quadratic equation
    a = F
    
    # Coefficient 'b' of the quadratic equation
    # This term represents the total initial energy of the system.
    b = 0.5 * M * v_0**2 - GmM / l_0 - F * l_0
    
    # Coefficient 'c' of the quadratic equation
    c = GmM

    # Calculate the discriminant (D = b^2 - 4ac) to find the roots.
    discriminant = b**2 - 4 * a * c

    # A real solution exists only if the discriminant is non-negative.
    # If it's negative, the spaceship has enough energy to escape to infinity.
    if discriminant < 0:
        print("The initial energy is too high for the spaceship to be trapped.")
        print("It will escape the asteroid's gravitational pull and travel to infinity.")
        final_answer = float('inf')
    else:
        # Calculate the two roots (potential turning points) of the quadratic equation
        sqrt_discriminant = math.sqrt(discriminant)
        root1 = (-b - sqrt_discriminant) / (2 * a)
        root2 = (-b + sqrt_discriminant) / (2 * a)

        # The spaceship starts at l_0 and moves away. The maximum distance is the
        # turning point it reaches first. This is the smaller root that is
        # greater than the initial distance l_0.
        l_max = min(root1, root2)
        
        print("To find the maximum distance (l_max), we solve the quadratic equation derived from the conservation of energy:")
        print(f"({a:.4e}) * l_max^2 + ({b:.4e}) * l_max + ({c:.4e}) = 0")
        print("\nThis equation gives two possible turning points for the spaceship:")
        print(f"Root 1: {max(root1, root2):.4e} m")
        print(f"Root 2: {min(root1, root2):.4e} m")
        print(f"\nSince the spaceship starts at l_0 = {l_0:.4e} m and moves outward, its motion is reversed at the first turning point it encounters.")
        print("\nTherefore, the maximum distance the spaceship can achieve is:")
        print(f"l_max = {l_max:.4e} meters")
        final_answer = l_max
    
    # Return final answer for submission format.
    return final_answer

if __name__ == '__main__':
    max_dist = calculate_max_distance()
    print(f"\n<<< {max_dist} >>>")
