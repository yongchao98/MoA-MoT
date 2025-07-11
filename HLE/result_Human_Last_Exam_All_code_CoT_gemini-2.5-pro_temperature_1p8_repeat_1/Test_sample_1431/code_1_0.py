import math

def calculate_lifting_force(m, l, v):
    """
    Calculates the force F required to lift a rope of mass m and length l
    such that its bottom end has speed v when it just leaves the ground.

    Args:
        m (float): Mass of the rope in kg.
        l (float): Length of the rope in meters.
        v (float): Final speed of the rope in m/s.
    """
    g = 9.8  # Acceleration due to gravity in m/s^2

    # The required force F is given by the formula: F = (m*g/2) + (m*v^2 / (2*l))
    # This is derived from the work-energy theorem: Work_done_by_F = Change_in_PE + Change_in_KE
    # F*l = m*g*(l/2) + 0.5*m*v^2

    potential_energy_term = (m * g) / 2
    kinetic_energy_term = (m * v**2) / (2 * l)
    
    force = potential_energy_term + kinetic_energy_term

    print("--- Mystical Force Calculation ---")
    print(f"Rope Mass (m): {m} kg")
    print(f"Rope Length (l): {l} m")
    print(f"Final Speed (v): {v} m/s")
    print(f"Gravity (g): {g} m/s^2")
    print("\nThe formula for the force F is: F = (m*g/2) + (m*v^2)/(2*l)")
    
    # Print the equation with the specific numbers
    print("\nSubstituting the values:")
    print(f"F = ({m} * {g} / 2) + ({m} * {v}**2 / (2 * {l}))")
    print(f"F = {potential_energy_term} + {kinetic_energy_term}")
    
    print(f"\nThe exact force F required is: {force:.2f} Newtons")

if __name__ == '__main__':
    # Example values for the monk's trial
    rope_mass = 15.0  # kg
    rope_length = 20.0  # meters
    final_speed = 4.0  # m/s
    
    calculate_lifting_force(rope_mass, rope_length, final_speed)