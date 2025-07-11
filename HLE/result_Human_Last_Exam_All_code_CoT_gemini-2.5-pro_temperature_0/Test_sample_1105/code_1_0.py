import math

def calculate_max_distance():
    """
    Calculates the maximum distance l_max a spaceship can achieve from an asteroid.
    
    The calculation is based on the work-energy theorem, resulting in a
    quadratic equation for l_max, which is then solved.
    """
    
    # --- Input Parameters ---
    # Mass of the asteroid A (kg)
    m = 5.972e12  # A small asteroid
    # Mass of the spaceship B (kg)
    M = 50000.0   # A large spaceship
    # Initial distance between A and B (m)
    l_0 = 2000.0
    # Initial speed of the spaceship (m/s)
    v_0 = 20.0
    # Additional constant force applied to the spaceship (N)
    F = 1000.0
    # Universal gravitational constant (m^3 kg^-1 s^-2)
    G = 6.67430e-11

    print("--- Input Parameters ---")
    print(f"Mass of asteroid (m): {m} kg")
    print(f"Mass of spaceship (M): {M} kg")
    print(f"Initial distance (l_0): {l_0} m")
    print(f"Initial speed (v_0): {v_0} m/s")
    print(f"Applied force (F): {F} N")
    print("-" * 26)

    # --- Calculation Steps ---
    
    # Step 1: Calculate the constant C from the work-energy equation.
    # C = F*l_0 + G*m*M/l_0 - 0.5*M*v_0^2
    term1 = F * l_0
    term2 = (G * m * M) / l_0
    term3 = 0.5 * M * v_0**2
    C = term1 + term2 - term3

    # Step 2: Calculate the term G*m*M, which is the constant term in the quadratic equation.
    GmM = G * m * M

    # Step 3: Calculate the discriminant of the quadratic equation.
    # Discriminant D = C^2 - 4*F*(G*m*M)
    try:
        discriminant = C**2 - 4 * F * GmM
    except ValueError:
        print("\nError: Calculation resulted in a math domain error, likely due to very large numbers.")
        return

    print("\n--- Solving the Quadratic Equation ---")
    print("The problem reduces to solving the quadratic equation for l_max:")
    print(f"a*(l_max)^2 + b*(l_max) + c = 0\n")
    print("Where the coefficients are:")
    print(f"a = F = {F}")
    print(f"b = -C = {-C}")
    print(f"c = G*m*M = {GmM}")
    
    # Step 4: Check if a real solution exists.
    if discriminant < 0:
        print("\nResult: The discriminant is negative.")
        print("This means the spaceship has enough energy to escape to infinity.")
        print("The maximum distance is infinite.")
    else:
        # Step 5: Calculate the two roots. The smaller root is the physical solution.
        sqrt_discriminant = math.sqrt(discriminant)
        
        # l_max = (C - sqrt(discriminant)) / (2*F)
        l_max = (C - sqrt_discriminant) / (2 * F)
        
        print(f"\n--- Final Result ---")
        print(f"The calculated maximum distance is:")
        print(f"l_max = {l_max:.4f} meters")

# Execute the calculation
calculate_max_distance()
<<<l_max = (C - math.sqrt(C**2 - 4*F*G*m*M)) / (2*F)>>>