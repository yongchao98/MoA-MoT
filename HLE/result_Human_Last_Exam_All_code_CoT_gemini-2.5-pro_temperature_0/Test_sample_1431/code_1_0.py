def solve_mystical_challenge(m, l, v, g=9.81):
    """
    Calculates the mystical force F required to lift a rope of mass m and length l,
    such that the top end has velocity v at the moment the rope is fully lifted.

    The formula used is F = m*g + (m*v^2)/(3*l).

    Args:
        m (float): mass of the rope in kilograms (kg).
        l (float): length of the rope in meters (m).
        v (float): final velocity of the top end in meters per second (m/s).
        g (float): acceleration due to gravity (m/s^2), default is 9.81.
    """
    if l <= 0:
        print("Error: The length of the rope 'l' must be positive.")
        return

    # Part 1: Force to support the rope's weight
    force_gravity = m * g

    # Part 2: Additional force due to the dynamics of the lift
    force_dynamic = (m * v**2) / (3 * l)

    # The total force is the sum of the two components
    total_force = force_gravity + force_dynamic

    print("The final equation for the mystical force F is: F = mg + (mv^2)/(3l)")
    print("\n--- Calculation Breakdown ---")
    print(f"Given parameters: m = {m} kg, l = {l} m, v = {v} m/s, g = {g} m/s^2")
    print("-" * 28)

    # Outputting each number in the final equation as requested
    print("1. Force to support the rope's weight (mg):")
    print(f"   F_gravity = {m} * {g}")
    print(f"   F_gravity = {force_gravity:.2f} N")
    print("\n2. Force required for the lifting dynamics (mv^2 / (3l)):")
    print(f"   F_dynamic = ({m} * {v}**2) / (3 * {l})")
    print(f"   F_dynamic = {force_dynamic:.2f} N")
    print("-" * 28)

    print("3. Total Mystical Force (F = F_gravity + F_dynamic):")
    print(f"   F = {force_gravity:.2f} N + {force_dynamic:.2f} N")
    print(f"   F = {total_force:.2f} N")
    print("\nThis is the exact force the monk must summon to succeed.")


# --- Example Challenge ---
# A monk must lift a rope with:
# Mass (m) = 15 kg
# Length (l) = 25 m
# Final speed (v) = 4 m/s

# Let's provide these values to our function
mass_rope = 15.0
length_rope = 25.0
final_speed = 4.0

solve_mystical_challenge(mass_rope, length_rope, final_speed)