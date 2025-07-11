import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    given an initial push, a constant applied force, and gravitational attraction.
    """
    # Define the constants and variables of the system.
    # G is the gravitational constant in m^3 kg^-1 s^-2
    G = 6.67408e-11 
    
    # Mass of the asteroid (m) in kg
    m = 2.158e17 
    # Mass of the spaceship (M) in kg
    M = 800.0
    # Initial distance (l_0) in meters
    l_0 = 8000.0
    # Initial velocity (v_0) in m/s
    v_0 = 2.0
    # Constant applied force (F) in Newtons
    F = 1.0

    print("Given physical parameters:")
    print(f"Mass of asteroid (m): {m:.2e} kg")
    print(f"Mass of spaceship (M): {M} kg")
    print(f"Initial distance (l_0): {l_0} m")
    print(f"Initial speed (v_0): {v_0} m/s")
    print(f"Applied force (F): {F} N")
    print(f"Gravitational Constant (G): {G:.3e} m^3 kg^-1 s^-2\n")

    # Step 1: Calculate the constant C from the initial conditions.
    # C = F*l_0 + G*M*m/l_0 - 0.5*M*v_0^2
    term_F_l0 = F * l_0
    term_G_Mm_l0 = (G * M * m) / l_0
    term_K_initial = 0.5 * M * v_0**2
    
    C = term_F_l0 + term_G_Mm_l0 - term_K_initial

    # Step 2: Form the quadratic equation: a*x^2 + b*x + c = 0
    # where x is l_max, a = F, b = -C, c = G*M*m
    a = F
    b = -C
    c = G * M * m
    
    print("The problem reduces to solving the quadratic equation for l_max:")
    # We output each number in the final equation as requested.
    print(f"{a} * l_max^2 + ({b:.4f}) * l_max + ({c:.4e}) = 0\n")


    # Step 3: Solve the quadratic equation.
    # First, calculate the discriminant (D = b^2 - 4ac)
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("The conditions are such that there is no real solution for the maximum distance.")
        print("This might happen if the initial energy is too low to move anywhere.")
        return

    # Calculate the two roots
    sqrt_discriminant = math.sqrt(discriminant)
    
    # Based on the physics of the problem, the smaller root corresponds to the maximum distance.
    # The spaceship starts at l_0, moves away, and stops at the first turning point it encounters.
    l_max_1 = (-b - sqrt_discriminant) / (2 * a)
    l_max_2 = (-b + sqrt_discriminant) / (2 * a)

    # The smaller root is the physical answer for l_max in this scenario.
    l_max = min(l_max_1, l_max_2)
    
    if l_max < l_0:
         print(f"Calculation indicates an issue: The calculated max distance ({l_max:.2f} m) is less than the initial distance ({l_0} m).")
         print("This setup might correspond to an escape trajectory under different initial velocity assumptions.")
    else:
        print(f"The maximum distance l_max is the smaller root of the quadratic equation.")
        print(f"Calculated Maximum Distance (l_max): {l_max:.2f} meters")

# Execute the function to solve the problem
calculate_max_distance()
<<<10000.00>>>