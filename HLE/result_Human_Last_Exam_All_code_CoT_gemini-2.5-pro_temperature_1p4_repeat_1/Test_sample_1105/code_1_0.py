import math

def calculate_max_distance(m, M, l0, v0, F, G=6.674e-11):
    """
    Calculates the maximum distance the spaceship can achieve from the asteroid.

    Args:
        m (float): Mass of the asteroid (kg)
        M (float): Mass of the spaceship (kg)
        l0 (float): Initial distance between A and B (m)
        v0 (float): Initial speed of the spaceship (m/s)
        F (float): Additional force applied to the spaceship (N)
        G (float): Gravitational constant (N m^2 / kg^2)
    """
    print("--- Input Parameters ---")
    print(f"Mass of asteroid (m): {m:.2e} kg")
    print(f"Mass of spaceship (M): {M:.2e} kg")
    print(f"Initial distance (l0): {l0:.2e} m")
    print(f"Initial speed (v0): {v0:.2e} m/s")
    print(f"Applied force (F): {F:.2e} N")
    print(f"Gravitational constant (G): {G:.2e} N m^2/kg^2")
    print("-" * 26)

    # Calculate G*M*m for convenience
    GMm = G * M * m

    # Condition for the ship to potentially slow down: net force must be negative at some point.
    # This occurs if l0 is less than the equilibrium distance where F_applied = F_gravity.
    l_eq_sq = GMm / F
    if l0**2 >= l_eq_sq:
        print("Condition for escape met: Initial force is already repulsive or zero.")
        print("The spaceship will continue to accelerate away and never stops.")
        print("Maximum distance is infinite.")
        return float('inf')

    # The equation for l_max is a quadratic: a*x^2 + b*x + c = 0
    # where x = l_max
    # a = F
    # b = -(F*l0 + GMm/l0 - 0.5*M*v0**2)
    # c = GMm
    
    # Let's calculate the term C = -b
    # C represents the value of the effective potential energy at the turning point.
    C = F * l0 + GMm / l0 - 0.5 * M * v0**2

    # For a real turning point to exist, C must be greater than the minimum
    # of the effective potential function H(l) = F*l + GMm/l, which is 2*sqrt(F*GMm).
    H_min = 2 * math.sqrt(F * GMm)

    if C < H_min:
        print(f"Condition for escape met: Initial kinetic energy is too high.")
        print(f"Required energy at turning point (C = {C:.4e}) is less than the minimum possible ({H_min:.4e}).")
        print("The spaceship's velocity never reaches zero.")
        print("Maximum distance is infinite.")
        return float('inf')

    # Calculate the discriminant D = b^2 - 4ac = C^2 - 4*F*GMm
    discriminant = C**2 - 4 * F * GMm

    # The two roots for l_max are (C +/- sqrt(D)) / (2*F)
    # The maximum distance reached (the first turning point) is the smaller root.
    l_max = (C - math.sqrt(discriminant)) / (2 * F)

    print("--- Calculation ---")
    print("The maximum distance l_max is the smaller root of the quadratic equation:")
    print("F * l_max^2 - C * l_max + GMm = 0")
    print("where C = F*l0 + G*M*m/l0 - 0.5*M*v0^2")
    print("\nCalculated values for the equation:")
    print(f"F = {F:.4e}")
    print(f"C = {C:.4e}")
    print(f"GMm = {GMm:.4e}")
    
    print("\nThe solution is l_max = (C - sqrt(C^2 - 4*F*GMm)) / (2*F)")
    print("\n--- Result ---")
    print(f"The maximum distance l_max is: {l_max:.4f} meters.")
    
    return l_max


# --- Example Calculation ---
# Let's define some sample values for the variables.
# These values are chosen such that a finite maximum distance exists.
m_asteroid = 1.0e18   # kg
M_spaceship = 1.0e6   # kg
l0_initial = 5.0e5    # m (500 km)
v0_initial = 1.0      # m/s
F_force = 100.0       # N

l_max_distance = calculate_max_distance(m_asteroid, M_spaceship, l0_initial, v0_initial, F_force)

# For the final answer format
if l_max_distance != float('inf'):
    final_answer = f"<<<{l_max_distance:.4f}>>>"
else:
    final_answer = "<<<inf>>>"

# The final line below would normally be hidden from the user, but we print it
# to show what would be returned in the requested format.
# print(final_answer)