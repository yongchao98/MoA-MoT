import math

def calculate_lifting_force(m, l, v):
    """
    Calculates the constant force F required to lift a rope of mass m and length l
    such that its final velocity is v when fully lifted.

    The formula is derived from the Work-Energy theorem:
    F = (m*g/2) + (m*v^2 / (2*l))
    """
    g = 9.8  # Acceleration due to gravity in m/s^2

    # --- Calculation Steps ---
    
    # 1. Calculate the component of the force needed to overcome gravity (change in potential energy)
    force_potential = (m * g) / 2
    
    # 2. Calculate the component of the force needed to achieve the final speed (change in kinetic energy)
    force_kinetic = (m * v**2) / (2 * l)
    
    # 3. The total force is the sum of the two components
    total_force = force_potential + force_kinetic

    # --- Output the results ---
    
    print("--- The Mystical Challenge of the Rope ---")
    print(f"Given values: mass (m) = {m} kg, length (l) = {l} m, final speed (v) = {v} m/s")
    print(f"Gravitational acceleration (g) is assumed to be {g} m/s^2.\n")
    
    print("The final equation for the force F is derived from the work-energy principle:")
    print("F * l = (m * g * l / 2) + (1/2 * m * v^2)")
    print("Solving for F: F = (m * g / 2) + (m * v^2 / (2 * l))\n")
    
    print("Let's substitute the numbers into the equation:")
    print(f"Force for potential energy: F_potential = ({m} * {g}) / 2 = {force_potential:.2f} N")
    print(f"Force for kinetic energy:  F_kinetic = ({m} * {v}^2) / (2 * {l}) = {force_kinetic:.2f} N\n")
    
    print("The total force F is the sum of these two parts:")
    print(f"F = F_potential + F_kinetic")
    print(f"F = {force_potential:.2f} + {force_kinetic:.2f}")
    
    final_force_str = f"{total_force:.2f}"
    print(f"The exact force F required is: {final_force_str} N")
    return final_force_str

# --- Example values for the challenge ---
# You can change these values to see how the force changes.
rope_mass_m = 10  # in kilograms
rope_length_l = 20 # in meters
final_speed_v = 5  # in meters/second

# Run the calculation and get the final answer
final_answer = calculate_lifting_force(rope_mass_m, rope_length_l, final_speed_v)

# The final answer in the required format
print(f"\n<<<F = {final_answer}>>>")