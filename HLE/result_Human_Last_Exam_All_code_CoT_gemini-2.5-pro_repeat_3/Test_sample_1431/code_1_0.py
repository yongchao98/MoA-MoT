import math

def calculate_lifting_force():
    """
    Calculates the force required to lift a rope of mass m and length l
    such that it has a speed v at the moment it's fully lifted.

    You can modify the values for m, l, and v below.
    """
    # --- Monk's Variables ---
    # Mass of the rope in kilograms (kg)
    m = 10.0
    # Length of the rope in meters (m)
    l = 20.0
    # Final speed of the rope in meters per second (m/s)
    v = 5.0

    # --- Physical Constant ---
    # Acceleration due to gravity in m/s^2
    g = 9.8

    # --- Calculation using the derived formula: F = (m/2) * (g + v^2/l) ---

    # Step 1: Calculate v^2 / l
    v_squared_over_l = v**2 / l

    # Step 2: Calculate g + v^2 / l
    parentheses_term = g + v_squared_over_l

    # Step 3: Calculate m / 2
    m_over_2 = m / 2

    # Step 4: Calculate the final force F
    force = m_over_2 * parentheses_term

    # --- Output the results step-by-step ---
    print("To solve the challenge, we use the work-energy theorem.")
    print("The formula for the force F is: F = (m/2) * (g + v^2/l)\n")

    print("Given values:")
    print(f"  Mass (m) = {m} kg")
    print(f"  Length (l) = {l} m")
    print(f"  Speed (v) = {v} m/s")
    print(f"  Gravity (g) = {g} m/s^2\n")

    print("Step-by-step calculation:")
    print(f"1. F = ({m} / 2) * ({g} + {v}^2 / {l})")
    print(f"2. F = ({m_over_2}) * ({g} + {v**2} / {l})")
    print(f"3. F = {m_over_2} * ({g} + {v_squared_over_l})")
    print(f"4. F = {m_over_2} * ({parentheses_term})")
    print(f"5. F = {force:.2f} N\n")

    print(f"The exact mystical force F required is {force:.2f} Newtons.")

# Execute the calculation
calculate_lifting_force()