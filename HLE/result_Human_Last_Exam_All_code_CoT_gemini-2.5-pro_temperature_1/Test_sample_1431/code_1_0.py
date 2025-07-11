import math

def calculate_lifting_force():
    """
    Calculates the force required to lift a rope of length l and mass m
    such that its end has speed v when it just leaves the ground.
    """
    # --- Parameters of the magical rope ---
    # You can change these values to match a specific problem.
    m = 10.0  # mass in kilograms (kg)
    l = 20.0  # length in meters (m)
    v = 5.0   # final speed in meters per second (m/s)
    
    # --- Physical Constant ---
    g = 9.8   # acceleration due to gravity in m/s^2

    print("--- The Monk's Strength Challenge ---")
    print("Calculating the mystical force F required to lift the rope.\n")

    print("Given Parameters:")
    print(f"Rope Mass (m): {m} kg")
    print(f"Rope Length (l): {l} m")
    print(f"Final Speed (v): {v} m/s")
    print(f"Gravity (g): {g} m/s^2\n")

    print("The solution is based on the Work-Energy Theorem: W_F = ΔK + ΔU")
    print("F * l = (1/2)*m*v^2 + m*g*(l/2)")
    print("Solving for F gives the formula:")
    print("F = (m * v^2) / (2 * l) + (m * g) / 2\n")

    # --- Calculation ---
    # Calculate the two terms in the equation for F
    kinetic_energy_term = (m * v**2) / (2 * l)
    potential_energy_term = (m * g) / 2

    # Calculate the total force F
    F = kinetic_energy_term + potential_energy_term

    print("--- Calculation with Numbers ---")
    # Output each number in the final equation
    print(f"F = ({m} * {v}**2) / (2 * {l}) + ({m} * {g}) / 2")
    print(f"F = {kinetic_energy_term:.2f} N (from kinetic energy) + {potential_energy_term:.2f} N (from potential energy)")
    print("-" * 30)
    print(f"Total Force (F) = {F:.2f} Newtons")
    print("-" * 30)

# Run the calculation
calculate_lifting_force()