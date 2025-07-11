import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    based on initial conditions and applied forces.
    """
    # --- Define physical constants and parameters ---
    # Gravitational constant in N m^2/kg^2
    G = 6.67430e-11
    # Mass of the asteroid (A) in kg (using a large value for a clear result)
    m = 6e24
    # Mass of the spaceship (B) in kg
    M = 5e4
    # Initial distance between A and B in meters
    l0 = 7e6
    # Initial speed of the spaceship in m/s
    v0 = 500
    # Additional constant force applied to the spaceship in Newtons
    F = 1e5

    # --- Solve for l_max using the work-energy theorem ---

    # The work-energy equation rearranges into a quadratic equation for l_max:
    # F * l_max^2 - (F*l0 + G*M*m/l0 - 0.5*M*v0^2) * l_max + G*M*m = 0
    # This is in the form a*x^2 + b*x + c = 0

    # Calculate the coefficients of the quadratic equation
    # a*x^2 + b*x + c = 0 where x = l_max
    a = F
    # The term 'C' from the derivation: F*l0 + G*M*m/l0 - 0.5*M*v0^2
    constant_C = (F * l0) + (G * M * m / l0) - (0.5 * M * v0**2)
    b = -constant_C
    c = G * M * m

    # Calculate the discriminant (d = b^2 - 4ac)
    discriminant = b**2 - 4 * a * c

    # Check if a solution exists
    if discriminant < 0:
        print("The spaceship has enough energy to escape to infinity.")
        print("There is no finite maximum distance.")
    else:
        # Calculate the two roots of the quadratic equation using the quadratic formula
        # x = (-b ± sqrt(discriminant)) / (2a)
        # The physically meaningful root for the maximum distance reached is the smaller one,
        # as it represents the first turning point.
        sqrt_discriminant = math.sqrt(discriminant)
        
        # In our derivation, b = -constant_C, so -b = constant_C
        l_max_1 = (constant_C - sqrt_discriminant) / (2 * a)
        l_max_2 = (constant_C + sqrt_discriminant) / (2 * a)
        
        # The smaller root is the correct one, provided the initial distance l0 is less than it.
        # Otherwise, the ship is already past the turning point and escapes.
        l_max = l_max_1

        if l0 >= l_max:
             print("The spaceship is already past its turning point and will escape to infinity.")
        else:
            print("The problem is solved using the work-energy theorem, which leads to the following quadratic equation for the maximum distance (l_max):")
            # Print the equation with the calculated coefficients
            print(f"({a:.4g}) * l_max^2 + ({b:.4g}) * l_max + ({c:.4g}) = 0")
            print("\nSolving this equation gives the maximum distance l_max.")
            print(f"The maximum distance l_max is: {l_max:,.2f} meters.")

# Execute the function
calculate_max_distance()