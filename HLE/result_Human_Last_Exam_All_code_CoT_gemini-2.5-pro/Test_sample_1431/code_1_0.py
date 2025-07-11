import math

def calculate_lifting_force():
    """
    Calculates the force required to lift a rope based on user-provided parameters.
    """
    # --- Introduction and Input ---
    print("--- Mystical Rope Force Calculator ---")
    print("This script calculates the force F required to lift a rope of mass m and length l,")
    print("such that its speed is v when it just leaves the ground.")
    print("\nPlease provide the following values:")
    
    try:
        m = float(input("Enter the mass of the rope (m) in kilograms: "))
        l = float(input("Enter the length of the rope (l) in meters: "))
        v = float(input("Enter the final speed of the rope (v) in meters/second: "))
    except ValueError:
        print("\nInvalid input. Please enter numeric values.")
        return

    if l <= 0:
        print("\nError: Length (l) must be a positive number.")
        return
        
    # --- Constants ---
    g = 9.8  # Acceleration due to gravity in m/s^2

    # --- Physics Explanation ---
    print("\n--- Calculation Breakdown ---")
    print("We use the work-energy theorem: Work_done = Change_in_Energy")
    print("Work done by force F is: W = F * l")
    print("Change in Potential Energy is: ΔPE = m * g * l / 2")
    print("Change in Kinetic Energy is: ΔKE = m * v^2 / 2")
    print("\nSetting W = ΔPE + ΔKE gives the formula:")
    print("F * l = (m * g * l / 2) + (m * v^2 / 2)")
    print("Solving for F, we get: F = (m * g / 2) + (m * v^2 / (2 * l))")

    # --- Calculation ---
    force_potential_part = (m * g) / 2
    force_kinetic_part = (m * v**2) / (2 * l)
    total_force = force_potential_part + force_kinetic_part

    # --- Output Final Equation with Numbers ---
    print("\nPlugging in your numbers:")
    print(f"m = {m} kg, l = {l} m, v = {v} m/s, g = {g} m/s^2")
    print("\nThe final equation is:")
    print(f"F = ({m} * {g} / 2) + ({m} * {v}**2 / (2 * {l}))")
    print(f"F = {force_potential_part:.2f} N (from potential energy) + {force_kinetic_part:.2f} N (from kinetic energy)")
    
    # --- Final Answer ---
    print("\n------------------------------------")
    print(f"The exact force F required is: {total_force:.2f} N")
    print("------------------------------------")


if __name__ == "__main__":
    calculate_lifting_force()
