import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    based on initial conditions and forces.
    """
    # --- Problem Parameters (using example values) ---
    # Gravitational constant
    G = 6.67430e-11  # m^3 kg^-1 s^-2
    # Mass of the asteroid (m)
    m = 5.972e20  # kg (a large asteroid)
    # Mass of the spaceship (M)
    M = 1.0e5  # kg (100 metric tons)
    # Initial distance (l_0)
    l_0 = 1.0e6  # meters (1000 km)
    # Initial speed (v_0)
    v_0 = 100.0  # m/s
    # Additional applied force (F)
    F = 1000.0  # Newtons

    print("--- Input Parameters ---")
    print(f"Gravitational Constant (G): {G:.4e} m^3 kg^-1 s^-2")
    print(f"Asteroid Mass (m): {m:.4e} kg")
    print(f"Spaceship Mass (M): {M:.4e} kg")
    print(f"Initial Distance (l_0): {l_0:.4e} m")
    print(f"Initial Speed (v_0): {v_0:.1f} m/s")
    print(f"Applied Force (F): {F:.1f} N\n")

    # --- Calculation Steps ---
    
    # 1. Calculate the constant C from the work-energy equation
    # C = F*l_0 + (G*m*M)/l_0 - (1/2)*M*v_0^2
    term1 = F * l_0
    term2 = (G * m * M) / l_0
    term3 = 0.5 * M * v_0**2
    C = term1 + term2 - term3
    
    # 2. Define coefficients for the quadratic equation: a*x^2 + b*x + c = 0
    # where x is l_max.
    a = F
    b = -C
    c = G * m * M

    print("--- Solving the Quadratic Equation for l_max ---")
    print("Equation form: a * l_max^2 + b * l_max + c = 0")
    print(f"a = F = {a}")
    print(f"b = - (F*l_0 + G*m*M/l_0 - 0.5*M*v_0^2) = {b:.4e}")
    print(f"c = G*m*M = {c:.4e}\n")

    # 3. Calculate the discriminant (delta = b^2 - 4ac)
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("The discriminant is negative. This means the spaceship has enough energy to escape to infinity.")
        print("There is no finite maximum distance.")
        l_max = float('inf')
    else:
        # 4. Use the quadratic formula to find the two roots for l_max.
        # The larger root is the maximum distance.
        sqrt_discriminant = math.sqrt(discriminant)
        l_max1 = (-b + sqrt_discriminant) / (2 * a)
        l_max2 = (-b - sqrt_discriminant) / (2 * a)
        
        # The physically relevant solution is the larger one, representing the outer turning point.
        l_max = max(l_max1, l_max2)

    print("--- Result ---")
    print(f"The maximum distance (l_max) the spaceship can achieve is: {l_max:.4e} meters.")
    
    # Returning the value for the final answer format
    return l_max

# Run the calculation and store the result
max_distance_result = calculate_max_distance()
print(f"\n<<< {max_distance_result} >>>")
