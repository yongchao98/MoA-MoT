import math

def calculate_charges(V, epsilon, L, a, b):
    """
    Calculates the total volume charge and surface charges on electrodes
    based on the formulas from the selected answer choice.

    Args:
        V (float): Applied DC voltage.
        epsilon (float): Constant permittivity of the material.
        L (float): Depth of the cylindrical electrodes.
        a (float): Radius of the inner electrode.
        b (float): Radius of the outer electrode.
    """
    if not (b > a > 0):
        print("Radii must be positive and b > a.")
        return

    # Common factor in the denominators
    denominator_factor = 1 - (a**2 / b**2)

    # Calculate total volume charge (q_v) based on Option B
    q_v_numerator = -4 * math.pi * L * epsilon * V
    q_v = q_v_numerator / denominator_factor

    # Calculate total surface charge on inner electrode (q_s_a) based on Option B
    q_s_a_numerator = 2 * math.pi * L * epsilon * V
    q_s_a = q_s_a_numerator / denominator_factor
    
    # Calculate total surface charge on outer electrode (q_s_b) based on Option B
    q_s_b_numerator = -4 * math.pi * L * epsilon * V * (a**2 / b**2)
    q_s_b = q_s_b_numerator / denominator_factor

    print("Based on the selected answer choice:")
    print(f"Total volume charge (q_v) = {q_v_numerator} / (1 - ({a}^2 / {b}^2)) = {q_v}")
    print(f"Total surface charge on inner electrode (q_s(r=a)) = {q_s_a_numerator} / (1 - ({a}^2 / {b}^2)) = {q_s_a}")
    print(f"Total surface charge on outer electrode (q_s(r=b)) = {q_s_b_numerator} / (1 - ({a}^2 / {b}^2)) = {q_s_b}")

# Example usage with symbolic values to show the structure
# Note: You can replace these with numerical values to get a numerical result.
V_val = 10.0  # Volts
epsilon_val = 8.854e-12  # F/m (permittivity of free space, as an example)
L_val = 0.5   # meters
a_val = 0.01  # meters
b_val = 0.02  # meters

calculate_charges(V_val, epsilon_val, L_val, a_val, b_val)