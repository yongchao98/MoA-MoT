import math

def calculate_max_distance():
    """
    Calculates the maximum distance the spaceship can achieve from the asteroid.
    """
    # Gravitational constant
    G = 6.67430e-11  # m^3 kg^-1 s^-2

    # --- User-defined parameters ---
    # Mass of the asteroid (m) in kg
    m = 1e15
    # Mass of the spaceship (M) in kg
    M = 5e4
    # Initial distance (l0) in meters
    l0 = 1000.0
    # Initial speed (v0) in m/s
    v0 = 150.0
    # Applied force (F) in Newtons
    F = 100000.0
    
    print("Given parameters:")
    print(f"  Mass of asteroid (m): {m:.2e} kg")
    print(f"  Mass of spaceship (M): {M:.2e} kg")
    print(f"  Initial distance (l0): {l0} m")
    print(f"  Initial speed (v0): {v0} m/s")
    print(f"  Applied force (F): {F} N")
    print("-" * 30)

    # Calculate the constant K from the work-energy equation
    # K = F*l0 + G*M*m/l0 - 0.5*M*v0^2
    term_F_l0 = F * l0
    term_G = (G * M * m) / l0
    term_KE = 0.5 * M * v0**2
    K = term_F_l0 + term_G - term_KE

    # For a finite maximum distance to exist, K must be positive,
    # and the discriminant of the resulting quadratic equation must be non-negative.
    
    if K <= 0:
        print("The spaceship has enough initial energy to escape.")
        print("Maximum distance is infinite.")
        return

    # The quadratic equation is: F * l_max^2 - K * l_max + (G*M*m) = 0
    # Let a = F, b = -K, c = G*M*m
    a = F
    b = -K
    c = G * M * m

    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("The spaceship's total energy is too high to be trapped.")
        print("Maximum distance is infinite.")
        return

    # Print the quadratic equation being solved: a*x^2 + b*x + c = 0
    print("The maximum distance (l_max) is a solution to the quadratic equation:")
    print(f"{a:.4e} * l_max^2 + ({b:.4e}) * l_max + {c:.4e} = 0")
    print("-" * 30)

    # The maximum distance is the smaller root of the quadratic equation.
    # l_max = (-b - sqrt(discriminant)) / (2*a)
    # which is equivalent to (K - sqrt(discriminant)) / (2*F)
    l_max = (K - math.sqrt(discriminant)) / (2 * F)

    print(f"The maximum distance the spaceship can achieve is:")
    print(f"l_max = {l_max:.4f} meters")


if __name__ == '__main__':
    calculate_max_distance()