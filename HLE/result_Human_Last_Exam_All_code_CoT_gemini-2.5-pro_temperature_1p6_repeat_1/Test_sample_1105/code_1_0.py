import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    based on initial conditions and forces.
    """
    # Define physical constants and parameters
    # Gravitational constant in m^3 kg^-1 s^-2
    G = 6.67430e-11
    # Mass of the asteroid (m) in kg
    m = 1.0e12
    # Mass of the spaceship (M) in kg
    M = 5.0e4
    # Initial distance (l_0) in meters
    l_0 = 1000.0
    # Initial speed (v_0) in m/s
    v_0 = 0.1
    # Applied constant force (F) in Newtons
    F = 1.0

    print("--- Given Parameters ---")
    print(f"Gravitational Constant (G): {G} m^3 kg^-1 s^-2")
    print(f"Asteroid Mass (m): {m:.1e} kg")
    print(f"Spaceship Mass (M): {M:.1e} kg")
    print(f"Initial Distance (l_0): {l_0} m")
    print(f"Initial Speed (v_0): {v_0} m/s")
    print(f"Applied Force (F): {F} N")
    print("-" * 26)

    # From the work-energy theorem, we derive the equation:
    # F*l_max + G*m*M/l_max = F*l_0 + G*m*M/l_0 - 0.5*M*v_0**2
    # Let K_c = F*l_0 + G*m*M/l_0 - 0.5*M*v_0**2
    # The equation becomes F*l_max^2 - K_c*l_max + G*m*M = 0

    # Calculate constant terms
    GmM = G * m * M
    initial_kinetic_energy = 0.5 * M * v_0**2

    # Calculate the constant K_c from the work-energy equation
    K_c = F * l_0 + GmM / l_0 - initial_kinetic_energy
    
    # The equation for the turning points l is: F*l^2 - K_c*l + GmM = 0
    # We solve this quadratic equation for l. The discriminant is D = K_c^2 - 4*F*GmM
    
    print("\n--- Calculation Steps ---")
    print("The final equation for l_max is derived from the work-energy theorem.")
    print("It takes the form of a quadratic equation: F * l_max^2 - K_c * l_max + (G*m*M) = 0")
    print(f"\nFirst, we calculate the constants for the equation:")
    print(f"Value of G*m*M = {G:.3e} * {m:.1e} * {M:.1e} = {GmM:.4g}")
    print(f"Value of K_c = ({F}*{l_0}) + ({GmM:.4g}/{l_0}) - (0.5*{M:.1e}*{v_0}^2) = {K_c:.4g}")
    
    # Calculate the discriminant
    discriminant = K_c**2 - 4 * F * GmM
    
    print(f"Discriminant (D) = K_c^2 - 4*F*(G*m*M) = {K_c:.4g}^2 - 4*{F}*{GmM:.4g} = {discriminant:.4g}")

    # Check if a real solution exists. If the discriminant is negative,
    # the spaceship has enough energy to escape to infinity.
    if discriminant < 0:
        print("\nThe discriminant is negative. The spaceship's velocity never reaches zero.")
        print("Result: The spaceship escapes to infinity, so there is no finite maximum distance.")
        return float('inf')
    else:
        # The solutions are given by the quadratic formula: l = (K_c ± sqrt(D)) / (2*F)
        # The maximum distance l_max is the larger root (the one with the '+').
        l_max = (K_c + math.sqrt(discriminant)) / (2 * F)

        print("\n--- Final Equation and Solution ---")
        print("Using the quadratic formula for the larger root:")
        print(f"l_max = (K_c + sqrt(D)) / (2 * F)")
        print(f"l_max = ({K_c:.4g} + sqrt({discriminant:.4g})) / (2 * {F})")
        print(f"l_max = ({K_c:.4g} + {math.sqrt(discriminant):.4g}) / {2 * F}")
        print(f"l_max = {K_c + math.sqrt(discriminant):.4g} / {2 * F}")
        print(f"\nResult: The maximum distance l_max is {l_max:.5g} meters.")
        return l_max

# Run the calculation and store the result
final_distance = calculate_max_distance()
print(f"\n<<<{final_distance:.5g}>>>")