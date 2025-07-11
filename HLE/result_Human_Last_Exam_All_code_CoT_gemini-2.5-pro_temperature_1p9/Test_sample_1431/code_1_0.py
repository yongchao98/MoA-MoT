def calculate_lifting_force(m, l, v):
    """
    Calculates the constant force F required to lift a rope of mass m and length l,
    such that its final speed is v when the end just leaves the ground.

    Args:
        m (float): Mass of the rope in kg.
        l (float): Length of the rope in meters.
        v (float): Final speed of the rope in m/s.
    """
    g = 9.8  # Gravitational acceleration in m/s^2

    print("--- The Mystical Challenge ---")
    print(f"Rope Mass (m): {m} kg")
    print(f"Rope Length (l): {l} m")
    print(f"Target Speed (v): {v} m/s")
    print(f"Gravitational Acceleration (g): {g} m/s^2")
    print("\n--- Applying the Work-Energy Theorem ---")

    # Formula: F = m*g/2 + m*v^2/(2*l)
    
    # Calculate final potential energy (PE)
    pe_final = m * g * (l / 2)
    print(f"1. Final Potential Energy (PE = m*g*l/2):")
    print(f"   PE = {m} * {g} * ({l} / 2) = {pe_final:.2f} Joules")

    # Calculate final kinetic energy (KE)
    ke_final = 0.5 * m * v**2
    print(f"2. Final Kinetic Energy (KE = 1/2*m*v^2):")
    print(f"   KE = 0.5 * {m} * {v}**2 = {ke_final:.2f} Joules")

    # Total work done by force F is Work = PE + KE
    total_work = pe_final + ke_final
    print(f"3. Total work done (W = PE + KE) must equal F * l:")
    print(f"   F * {l} = {pe_final:.2f} + {ke_final:.2f} = {total_work:.2f} Joules")

    # Calculate the force F
    force = total_work / l
    
    print("\n--- The Required Force (F) ---")
    print(f"The formula for the force F is: F = (m/2) * (g + v^2/l)")
    print(f"Plugging in the values:")
    # The final equation is F = mg/2 + mv^2/(2l)
    print(f"F = ({m} * {g} / 2) + ({m} * {v}**2 / (2 * {l}))")
    term1 = m*g/2
    term2 = m*v**2/(2*l)
    print(f"F = {term1:.2f} + {term2:.2f}")
    print(f"F = {force:.2f} Newtons")

# --- Monk's Inputs ---
# You can change these values to test with a different rope
rope_mass = 12.0  # in kg
rope_length = 10.0 # in meters
target_speed = 5.0  # in m/s

# Run the calculation
calculate_lifting_force(rope_mass, rope_length, target_speed)

# To find the answer, let's derive it algebraically from the equation F = mg/2 + mv^2/(2l)
# This can be factored into F = (m/2) * (g + v^2/l)
final_force = (rope_mass/2)*(9.8 + target_speed**2 / rope_length)
<<<F = (m/2) * (g + v^2/l)>>>