def calculate_mystical_force(m, l, v):
    """
    Calculates the force required to lift a rope of mass 'm' and length 'l'
    such that the top end has speed 'v' at the moment the bottom end leaves the ground.

    Args:
        m (float): The mass of the rope in kilograms.
        l (float): The length of the rope in meters.
        v (float): The speed of the lifted end in meters per second.
    
    Returns:
        float: The required force in Newtons.
    """
    # Gravitational acceleration constant
    g = 9.8  # m/s^2

    # The force F has two components:
    # 1. The force to support the rope's weight (mg)
    # 2. The force required to impart momentum to the rope sections as they are lifted
    # The final formula is F = m*g + (m*v^2)/(2*l)

    weight_component = m * g
    momentum_component = (m * v**2) / (2 * l)
    total_force = weight_component + momentum_component

    # Printing the detailed calculation as requested
    print("The mystical force F is the sum of the rope's weight and a dynamic component related to its motion.")
    print(f"Formula: F = (m * g) + (m * v^2) / (2 * l)")
    print(f"Substituting the values m={m} kg, l={l} m, v={v} m/s, g={g} m/s^2:")
    print(f"F = ({m} * {g}) + ({m} * {v}^2) / (2 * {l})")
    print(f"F = {weight_component:.2f} + {momentum_component:.2f}")
    print(f"F = {total_force:.2f} Newtons")

    return total_force

if __name__ == '__main__':
    # Example values for the magical rope challenge
    rope_mass = 10.0  # kilograms
    rope_length = 20.0  # meters
    final_speed = 5.0  # meters per second

    calculate_mystical_force(rope_mass, rope_length, final_speed)
