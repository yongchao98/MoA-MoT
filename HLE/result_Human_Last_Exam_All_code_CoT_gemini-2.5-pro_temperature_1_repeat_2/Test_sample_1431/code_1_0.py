import sympy

def solve_rope_challenge():
    """
    Calculates the force F required to lift a rope of mass m and length l,
    such that the top end has speed v when the bottom end leaves the ground.
    
    The formula derived from physics principles is:
    F = m*g + (m*v^2) / (2*l)
    
    This script will print the formula and then compute an example.
    """
    
    # Define symbolic variables for the formula
    m, g, v, l = sympy.symbols('m g v l')
    
    # The formula for the force F
    force_formula = m*g + (m*v**2)/(2*l)
    
    print("The challenge is to find the force F. Based on the principles of mechanics for a variable-mass system, the formula is:")
    print(f"F = {force_formula}")
    print("\n-----------------------------------\n")

    # Let's use some example values for the monk's trial.
    # We will use g = 9.8 m/s^2 (standard gravity)
    example_m = 10.0  # kg
    example_l = 20.0  # meters
    example_v = 5.0   # m/s
    example_g = 9.8   # m/s^2

    # Calculate the components of the force
    weight_component = example_m * example_g
    kinetic_component = (example_m * example_v**2) / (2 * example_l)
    
    # Calculate the total force
    total_force = weight_component + kinetic_component
    
    print("To demonstrate, let's use the following example values:")
    print(f"Mass (m) = {example_m} kg")
    print(f"Length (l) = {example_l} m")
    print(f"Speed (v) = {example_v} m/s")
    print(f"Gravity (g) = {example_g} m/s^2")
    print("\nNow, let's calculate the final force by substituting these values into the formula:")
    
    # Print the equation with the numbers plugged in
    print(f"\nFinal Equation:")
    print(f"F = ({example_m} * {example_g}) + ({example_m} * {example_v}**2) / (2 * {example_l})")
    print(f"F = {weight_component} + {kinetic_component}")
    print(f"F = {total_force} Newtons")

# Execute the function to solve the challenge
solve_rope_challenge()
