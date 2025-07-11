import math

def calculate_max_distance(m, M, l0, v0, F, G=6.67430e-11):
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid.

    Args:
        m (float): Mass of the asteroid (kg)
        M (float): Mass of the spaceship (kg)
        l0 (float): Initial distance between A and B (m)
        v0 (float): Initial speed of the spaceship (m/s)
        F (float): Additional constant force on the spaceship (N)
        G (float): Gravitational constant (m^3 kg^-1 s^-2)
    """
    print(f"Calculating for parameters:")
    print(f"  m = {m} kg, M = {M} kg, l0 = {l0} m, v0 = {v0} m/s, F = {F} N, G = {G}")
    print("-" * 30)

    # From the work-energy theorem, we derive the equation:
    # F * l^2 - K_prime * l + G*m*M = 0
    # where K_prime = F*l0 + G*m*M/l0 - 0.5*M*v0^2
    
    # Calculate key terms
    gmm = G * m * M
    k_prime = F * l0 + gmm / l0 - 0.5 * M * v0**2
    g_min = 2 * math.sqrt(F * gmm)

    # Check for escape condition. If K' is less than the minimum of g(l),
    # the spaceship's kinetic energy never reaches zero.
    if k_prime < g_min:
        print("The spaceship has enough energy to escape to infinity.")
        print(f"This is because the energy term K' ({k_prime:.4f}) is less than the required minimum barrier ({g_min:.4f}).")
        print("Final Equation: No real solution for l_max.")
        print("l_max = Infinity")
        return

    # If we are here, a finite l_max exists. Solve the quadratic equation.
    # A*x^2 + B*x + C = 0 where A=F, B=-k_prime, C=gmm
    discriminant = k_prime**2 - 4 * F * gmm

    # The maximum distance is the smaller of the two positive roots.
    l_max = (k_prime - math.sqrt(discriminant)) / (2 * F)

    # Output the equation with the numbers filled in.
    print("The maximum distance (l_max) is the smaller positive root of the quadratic equation:")
    print(f"Final Equation: ({F}) * l_max^2 - ({k_prime:.4f}) * l_max + ({gmm:.4f}) = 0")
    print(f"\nThe calculated maximum distance is: {l_max:.4f} meters.")


if __name__ == "__main__":
    # --- Example 1: Working case with a finite distance ---
    # These parameters are chosen to ensure a finite solution exists.
    print("--- Example 1: Finite Distance ---")
    calculate_max_distance(m=1, M=1, l0=0.5, v0=0.5, F=1, G=1)
    
    print("\n" + "="*50 + "\n")

    # --- Example 2: Escape case with realistic numbers ---
    # Parameters for a more realistic scenario
    print("--- Example 2: Escape to Infinity ---")
    calculate_max_distance(
        m=6.687e15,  # Mass of a small asteroid (kg)
        M=5000,      # Mass of a spaceship (kg)
        l0=100000,   # Initial distance of 100 km (m)
        v0=10,       # Initial speed of 10 m/s
        F=0.2        # A small applied force (N)
    )