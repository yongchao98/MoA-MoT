import math

def calculate_charges(V, epsilon, L, a, b):
    """
    Calculates the total volume and surface charges based on the formulas in Option B.

    Args:
        V (float): DC voltage.
        epsilon (float): Permittivity of the material.
        L (float): Depth of the electrodes.
        a (float): Radius of the inner electrode.
        b (float): Radius of the outer electrode.
    """
    pi = math.pi
    
    # Check for valid radii
    if a >= b:
        print("Error: Radius 'a' must be less than radius 'b'.")
        return

    # Common denominator term
    denom_factor = (1 - (a**2 / b**2))

    # Formulas from Option B
    q_v_B = (-4 * V * epsilon * pi * L) / denom_factor
    q_s_a_B = (2 * pi * L * V * epsilon) / denom_factor
    q_s_b_B = (-4 * pi * L * V * epsilon * a**2) / (b**2 * denom_factor)

    print("Based on the formulas in Option B:")
    print(f"Total volume charge (q_v) = {q_v_B}")
    print(f"Total surface charge on inner electrode (q_s(r=a)) = {q_s_a_B}")
    print(f"Total surface charge on outer electrode (q_s(r=b)) = {q_s_b_B}")
    
    # For demonstration, let's check charge neutrality for Option B's formulas
    total_charge_B = q_v_B + q_s_a_B + q_s_b_B
    print(f"\nCheck for charge neutrality (q_v + q_s(a) + q_s(b)): {total_charge_B}")
    if not math.isclose(total_charge_B, 0):
        print("Note: The formulas in Option B are not self-consistent and do not satisfy charge neutrality (sum should be 0).")

# Example usage with placeholder values
# Let's use SI units for clarity
V_val = 100.0  # Volts
epsilon_val = 8.854e-12 * 2.5 # Permittivity of a dielectric (F/m)
L_val = 0.5    # meters
a_val = 0.01   # meters (1 cm)
b_val = 0.02   # meters (2 cm)

calculate_charges(V_val, epsilon_val, L_val, a_val, b_val)
