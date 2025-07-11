import math

def calculate_lifting_force():
    """
    Calculates the constant force F required to lift a rope of mass m and length l,
    such that its final speed is v when its end just leaves the ground.
    The function also prints the step-by-step calculation.
    """
    # Parameters of the mystical rope and challenge (example values)
    m = 10.0  # mass in kg
    l = 20.0  # length in meters
    v = 5.0   # final speed in m/s
    g = 9.8   # acceleration due to gravity in m/s^2

    # Explain the given values
    print("--- Solving the Mystical Challenge ---")
    print(f"Assuming the following values:")
    print(f"Mass of the rope (m): {m} kg")
    print(f"Length of the rope (l): {l} m")
    print(f"Final speed (v): {v} m/s")
    print(f"Gravity (g): {g} m/s^2")
    print("-" * 35)

    # --- Calculation based on the Work-Energy Theorem ---
    # The formula for the constant force F is: F = (m*g)/2 + (m*v**2)/(2*l)
    print("The formula for the required constant force F is:")
    print("F = (m * g / 2) + (m * v^2) / (2 * l)")
    print("")

    # Calculate the two terms separately for clarity
    # Term 1 corresponds to the work needed to overcome gravity (increase potential energy)
    # Term 2 corresponds to the work needed to give the rope its final speed (increase kinetic energy)
    potential_term = (m * g) / 2
    kinetic_term = (m * v**2) / (2 * l)

    # Print the equation with numbers plugged in
    print("Substituting the values into the equation:")
    print(f"F = ({m} * {g} / 2) + ({m} * {v}**2) / (2 * {l})")
    
    # Print the value of each term in the sum
    print(f"F = {potential_term} N + {kinetic_term} N")

    # Calculate the final force F
    F = potential_term + kinetic_term

    # Print the final result
    print("-" * 35)
    print(f"The exact force F required is: {F} Newtons.")
    print("-" * 35)


# Execute the function to solve the challenge
calculate_lifting_force()