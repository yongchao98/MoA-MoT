import math

def calculate_monk_force(m, l, v):
    """
    Calculates the mystical force F required to lift a rope of mass m and length l.

    Args:
        m (float): The mass of the rope in kilograms.
        l (float): The length of the rope in meters.
        v (float): The final speed of the rope's end in meters per second.
    """
    # Standard acceleration due to gravity in m/s^2
    g = 9.8

    if m <= 0 or l <= 0:
        print("Mass and length must be positive values.")
        return

    # The force F is the sum of the force to counteract gravity (mg)
    # and the force to give momentum to the rope as it's lifted.
    # The derived formula is F = m*g + (m*v^2)/(2*l).

    force_gravity = m * g
    force_momentum = (m * v**2) / (2 * l)
    total_force = force_gravity + force_momentum

    print("The mystical challenge requires calculating the force F using the formula:")
    print("F = m*g + (m*v^2)/(2*l)")
    print("\nGiven values:")
    print(f"  Mass (m): {m} kg")
    print(f"  Length (l): {l} m")
    print(f"  Speed (v): {v} m/s")
    print(f"  Gravity (g): {g} m/s^2")
    print("\nCalculation steps:")
    print(f"  Force against gravity = {m} * {g} = {force_gravity:.2f} N")
    print(f"  Force for momentum = ({m} * {v}**2) / (2 * {l}) = {force_momentum:.2f} N")
    print("\nFinal Equation:")
    # Using the print function to explicitly show the numbers in the equation
    print(f"F = ({m} * {g}) + ({m} * {v}**2) / (2 * {l})")
    print(f"F = {force_gravity:.2f} + {force_momentum:.2f}")
    print(f"F = {total_force:.2f} N")


# --- Example Usage ---
# You can change these values to test with different ropes.
rope_mass = 15.0  # kg
rope_length = 30.0  # meters
final_speed = 4.0   # m/s

calculate_monk_force(rope_mass, rope_length, final_speed)