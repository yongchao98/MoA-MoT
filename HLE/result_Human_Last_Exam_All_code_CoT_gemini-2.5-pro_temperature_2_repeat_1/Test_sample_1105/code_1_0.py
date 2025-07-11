import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    based on initial conditions and forces.
    """
    # Define constants and initial parameters.
    # These values are chosen to ensure a physically realistic scenario with a finite max distance.
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    m = 1.0e18       # Mass of the asteroid (m) in kg
    M = 1000.0       # Mass of the spaceship (M) in kg
    l_0 = 1000.0     # Initial distance (l_0) in meters
    v_0 = 100.0      # Initial speed (v_0) in m/s
    F = 10.0         # Applied force (F) in Newtons

    print("Solving for l_max based on the quadratic equation: a*(l_max)^2 + b*l_max + c = 0\n")

    # Calculate the coefficients of the quadratic equation
    # a = F
    # b = 1/2*M*v_0^2 - F*l_0 - G*m*M/l_0
    # c = G*m*M

    a = F
    ke_initial_term = 0.5 * M * v_0**2
    force_potential_term = F * l_0
    gravity_potential_term = (G * m * M) / l_0
    
    b = ke_initial_term - force_potential_term - gravity_potential_term
    c = G * m * M

    print(f"--- Equation Coefficients ---")
    print(f"Coefficient a (from F): {a}")
    print(f"Coefficient b (from initial total energy): {b:.4g}")
    print(f"Coefficient c (from G*m*M): {c:.4g}")
    print("-" * 27 + "\n")

    # Calculate the discriminant (b^2 - 4ac)
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("The spaceship has enough energy to escape to infinity.")
        print("There is no finite maximum distance.")
    else:
        # Calculate the two roots using the quadratic formula.
        # The larger root is the maximum distance l_max.
        sqrt_discriminant = math.sqrt(discriminant)
        l_max = (-b + sqrt_discriminant) / (2 * a)
        
        # We also calculate the other root for completeness, which is the inner turning point.
        # l_inner = (-b - sqrt_discriminant) / (2 * a)

        print("--- Final Calculation ---")
        # To satisfy the "output each number in the final equation" requirement,
        # we will show the values used in the quadratic formula calculation.
        print(f"l_max = ( -({b:.4g}) + sqrt( ({b:.4g})^2 - 4*({a})*({c:.4g}) ) ) / ( 2*({a}) )")
        print(f"\nThe calculated maximum distance (l_max) is: {l_max:.2f} meters")

# Execute the calculation
calculate_max_distance()