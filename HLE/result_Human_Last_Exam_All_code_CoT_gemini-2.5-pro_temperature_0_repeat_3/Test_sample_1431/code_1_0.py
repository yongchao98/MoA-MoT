import math

def calculate_mystical_force(m, l, v):
    """
    Calculates the force F required to lift a rope of mass m and length l,
    such that its speed is v at the moment it fully leaves the ground.

    Args:
        m (float): Mass of the rope in kg.
        l (float): Length of the rope in meters.
        v (float): Final speed of the rope in m/s.
    """
    # Standard acceleration due to gravity in m/s^2
    g = 9.8

    # The force F is the sum of two components:
    # 1. The force to counteract gravity for the entire rope (mg).
    # 2. The force to provide momentum to the part of the rope being lifted (mv^2/l).
    
    force_gravity = m * g
    force_impulse = (m * v**2) / l
    
    total_force = force_gravity + force_impulse
    
    # Print the explanation and the final equation with all numbers
    print("The total force F is the sum of the gravitational force and the impulse force.")
    print(f"F = F_gravity + F_impulse")
    print(f"F = (m * g) + (m * v^2 / l)")
    print(f"F = ({m} kg * {g} m/s^2) + ({m} kg * {v}^2 m/s^2 / {l} m)")
    print(f"F = {force_gravity:.2f} N + {force_impulse:.2f} N")
    print(f"The exact force F required is: {total_force:.2f} N")

if __name__ == '__main__':
    # You can change these values to solve for different ropes and speeds.
    # Example values for the mystical rope:
    rope_mass = 15.0  # kg
    rope_length = 20.0  # meters
    final_speed = 4.0  # m/s
    
    calculate_mystical_force(rope_mass, rope_length, final_speed)