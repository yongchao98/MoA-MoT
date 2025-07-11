import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    based on initial conditions and applied forces.
    """
    # Gravitational constant (in m^3 kg^-1 s^-2)
    G = 6.67430e-11

    # --- Input Parameters ---
    # To ensure a valid scenario where a maximum distance exists, we'll use a
    # consistent set of parameters. The spaceship starts near an "inner" turning
    # point and is given a push to travel to an "outer" turning point.

    # Mass of the asteroid (kg)
    m = 6e24 # e.g., Earth-like mass
    # Mass of the spaceship (kg)
    M = 5e4
    # Initial distance between A and B (m)
    l0 = 6.5e6
    # Initial speed of B away from A (m/s)
    v0 = 1000
    # Additional constant force on B (N)
    F = 1e5

    print(f"Given Parameters:")
    print(f"  Asteroid Mass (m): {m:.2e} kg")
    print(f"  Spaceship Mass (M): {M:.2e} kg")
    print(f"  Initial Distance (l0): {l0:.2e} m")
    print(f"  Initial Speed (v0): {v0:.2e} m/s")
    print(f"  Applied Force (F): {F:.2e} N")
    print("-" * 30)

    # Calculate initial kinetic energy
    K_initial = 0.5 * M * v0**2

    # We derived the quadratic equation: a*l_max^2 + b*l_max + c = 0
    # where a = F, c = G*M*m, and b = -C
    # C = F*l0 + (G*M*m)/l0 - K_initial

    GmM = G * M * m

    # Calculate the coefficients of the quadratic equation
    a = F
    C = F * l0 + GmM / l0 - K_initial
    b = -C
    c = GmM

    print("The problem is solved using the quadratic equation: a*x^2 + b*x + c = 0, where x is l_max.")
    print(f"The calculated coefficients are:")
    print(f"  a = {a}")
    print(f"  b = {b}")
    print(f"  c = {c}")
    print("-" * 30)


    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        # This case corresponds to an escape trajectory where the spaceship
        # has too much energy to ever turn around.
        print("The spaceship has escape velocity and will travel to infinity.")
        print("No maximum distance exists.")
        return

    # Solve for the two roots (turning points) using the quadratic formula
    sqrt_discriminant = math.sqrt(discriminant)
    
    # The larger root corresponds to the maximum distance l_max
    l_max = (-b + sqrt_discriminant) / (2 * a)

    print(f"The maximum distance the spaceship can achieve is:")
    print(f"l_max = {l_max:.4f} meters")
    print(f"l_max = {l_max/1000:.4f} kilometers")


# Run the calculation
calculate_max_distance()