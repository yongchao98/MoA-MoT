import math

def calculate_lifting_force(m, l, v):
    """
    Calculates the force F required to lift a rope of mass m and length l,
    such that its end has speed v the moment it leaves the ground.

    Args:
        m (float): Mass of the rope in kg.
        l (float): Length of the rope in meters.
        v (float): Final speed of the rope's end in m/s.
    """
    # Standard gravity constant in m/s^2
    g = 9.8

    # The force F is composed of two parts:
    # 1. The force to counteract the rope's full weight.
    weight_component = m * g
    
    # 2. The force required to accelerate new sections of the rope from rest to speed v.
    # This is derived from the variable-mass motion equation (F = v * dm/dt).
    kinetic_component = (m * v**2) / l

    # The total force is the sum of these two components.
    # This assumes the acceleration of the rope is zero at the final moment.
    total_force = weight_component + kinetic_component

    print("The challenge requires calculating the force F at the precise moment the rope leaves the ground.")
    print("The formula for this force is: F = m*g + (m*v^2)/l\n")
    print("Given values:")
    print(f"  Mass (m): {m} kg")
    print(f"  Length (l): {l} m")
    print(f"  Final Speed (v): {v} m/s")
    print(f"  Gravity (g): {g} m/s^2\n")

    print("Step-by-step calculation:")
    # Using 'int' or 'float' to display cleaner numbers if they don't have decimals.
    m_str = f"{int(m)}" if m == int(m) else f"{m}"
    l_str = f"{int(l)}" if l == int(l) else f"{l}"
    v_str = f"{int(v)}" if v == int(v) else f"{v}"
    
    # Print the full equation with the numbers plugged in
    print(f"F = ({m_str} * {g}) + ({m_str} * {v_str}^2) / {l_str}")
    print(f"F = {weight_component:.2f} + {kinetic_component:.2f}")
    print(f"F = {total_force:.2f} N")
    
    print("\nTherefore, the mystical force F required is {:.2f} Newtons.".format(total_force))


# --- Example values for the mystical rope ---
# You can change these values to see the result for a different rope.
rope_mass = 15.0      # in kilograms
rope_length = 25.0    # in meters
final_speed = 4.0     # in meters per second

calculate_lifting_force(rope_mass, rope_length, final_speed)
