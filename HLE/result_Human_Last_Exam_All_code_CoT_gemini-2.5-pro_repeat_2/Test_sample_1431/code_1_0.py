import math

def calculate_lifting_force():
    """
    Calculates the force required to lift a rope of length l and mass m
    such that it reaches a speed v just as it leaves the ground.
    
    The user can modify the example values below.
    """
    
    # --- Example Values (you can change these) ---
    l = 15.0  # length of the rope in meters
    m = 20.0  # mass of the rope in kilograms
    v = 5.0   # final speed in m/s
    # ---------------------------------------------

    # Gravitational constant
    g = 9.8  # m/s^2

    # Using the work-energy theorem: F * l = (1/2)*m*v^2 + m*g*(l/2)
    # Solving for F: F = (m*v^2)/(2*l) + (m*g)/2
    force = (m * v**2) / (2 * l) + (m * g) / 2

    # --- Output ---
    print("--- Calculating the Mystical Force F ---")
    print(f"Given values: rope length (l) = {l} m, mass (m) = {m} kg, final speed (v) = {v} m/s.")
    print(f"We use the formula: F = (m * v^2) / (2 * l) + (m * g) / 2\n")
    
    # Show the equation with numbers substituted
    print("Final Equation with numbers:")
    print(f"F = ({m} * {v}**2) / (2 * {l}) + ({m} * {g}) / 2")
    
    # Print the calculated result
    print(f"\nThe calculated force F required is: {force:.2f} Newtons")

# Execute the function
calculate_lifting_force()