import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    given an initial velocity, a constant applied force, and gravitational pull.
    """
    # Define the physical constants and parameters with placeholder values.
    # You can change these values to match a specific scenario.
    
    # Gravitational constant in N*(m/kg)^2
    G = 6.67430e-11
    # Mass of the asteroid A in kg
    m = 1.0e12
    # Mass of the spaceship B in kg
    M = 5.0e5
    # Initial distance between A and B in meters
    l0 = 1000.0
    # Initial speed of the spaceship in m/s
    v0 = 50.0
    # Additional constant force applied to the spaceship in Newtons
    F = 10000.0

    # The problem is solved by finding the roots of a quadratic equation derived
    # from the work-energy theorem: a*x^2 + b*x + c = 0, where x = l_max.
    # a = F
    # b = -(F*l0 + (G*m*M)/l0 - 0.5*M*v0^2)
    # c = G*m*M
    
    # Calculate the coefficients of the quadratic equation
    a = F
    GmM = G * m * M
    C = F * l0 + GmM / l0 - 0.5 * M * v0**2
    b = -C
    c = GmM
    
    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c
    
    print("This problem is solved by finding the maximum distance 'l_max' that satisfies the energy conservation equation.")
    print("\nThe equation is of the form: F * l_max^2 - C * l_max + G*m*M = 0")
    print("\nSubstituting the given values, the equation becomes:")
    print(f"({a:.1f}) * l_max^2 + ({b:.1f}) * l_max + ({c:.1f}) = 0\n")

    # Check if a real solution exists
    if discriminant < 0:
        print("The spaceship has enough energy to escape to infinity.")
        print("There is no finite maximum distance.")
        l_max = float('inf')
    else:
        # The quadratic formula gives two solutions for the distance.
        # The larger root is the maximum distance l_max.
        sol1 = (-b - math.sqrt(discriminant)) / (2 * a)
        sol2 = (-b + math.sqrt(discriminant)) / (2 * a)
        
        l_max = max(sol1, sol2)
        
        print(f"The calculated maximum distance l_max is: {l_max:.2f} meters")

    return l_max

if __name__ == '__main__':
    max_distance = calculate_max_distance()
    # The final answer is the numerical value calculated above.
    print(f"\n<<< {max_distance:.2f} >>>")
