import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    based on initial conditions and applied forces.
    """
    # --- User Inputs ---
    # Please replace these placeholder values with the actual values from your problem.
    m = 1.0e15  # mass of the asteroid A in kg
    M = 5.0e5   # mass of the spaceship B in kg
    l0 = 1.0e5  # initial distance between A and B in meters
    v0 = 100.0  # initial speed of the spaceship in m/s
    F = 100.0   # additional applied force in Newtons

    # Gravitational constant
    G = 6.67430e-11  # m^3 kg^-1 s^-2

    # The problem assumes that a finite maximum distance exists.
    # This is true if the initial push isn't strong enough to escape completely.
    # This requires l0 < sqrt(G*M*m/F) and v0 to be below the escape threshold.

    # The work-energy theorem leads to the quadratic equation:
    # F * l_max^2 - K * l_max + G*M*m = 0
    # where K is a constant derived from the initial conditions.

    # Calculate the intermediate constant K
    try:
        K = F * l0 + (G * M * m) / l0 - 0.5 * M * v0**2
    except ZeroDivisionError:
        print("Error: Initial distance l0 cannot be zero.")
        return

    # Calculate the terms of the quadratic equation a*x^2 + b*x + c = 0
    a = F
    b = -K
    c = G * M * m

    # Calculate the discriminant of the quadratic equation
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("The spaceship has enough energy to escape to infinity.")
        print("A finite maximum distance does not exist under these conditions.")
        return

    # The two roots of the quadratic equation represent the turning points.
    # The physically correct solution for the maximum distance reached is the smaller root.
    l_max = (-b - math.sqrt(discriminant)) / (2 * a)

    # --- Output the results ---
    print("--- Calculation of Maximum Distance (l_max) ---")
    print(f"Based on the work-energy principle, we solve for l_max in the equation:")
    print(f"F*l_max^2 - K*l_max + G*M*m = 0\n")
    
    print("Given Parameters:")
    print(f"  Mass of asteroid (m): {m:.2e} kg")
    print(f"  Mass of spaceship (M): {M:.2e} kg")
    print(f"  Initial distance (l0): {l0:.2e} m")
    print(f"  Initial speed (v0): {v0:.2f} m/s")
    print(f"  Applied force (F): {F:.2f} N")
    print(f"  Gravitational Constant (G): {G:.5e} m^3 kg^-1 s^-2\n")

    print("The final calculation is:")
    # Using f-strings to embed the numbers directly into the equation string
    print(f"l_max = (K - sqrt(K^2 - 4*F*G*M*m)) / (2*F)")
    print(f"l_max = ({K:.2e} - sqrt({K:.2e}^2 - 4 * {F:.2f} * {G:.5e} * {M:.2e} * {m:.2e})) / (2 * {F:.2f})")
    
    print("\n--- Result ---")
    print(f"The maximum distance l_max is: {l_max:.4e} meters.")
    
    # Return the final value for verification
    return l_max

# Run the calculation and store the result
max_distance = calculate_max_distance()
# The final answer is printed by the function.
# To comply with the output format, we print the numerical result here again.
if max_distance is not None:
    print(f"\n<<< {max_distance} >>>")
