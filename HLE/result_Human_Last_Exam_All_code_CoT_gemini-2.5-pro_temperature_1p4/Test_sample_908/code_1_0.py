import math

def calculate_charges(V, epsilon, L, a, b):
    """
    Calculates the total volume and surface charges based on the formulas in Option B.

    Args:
        V (float): Applied DC voltage.
        epsilon (float): Permittivity of the material.
        L (float): Length of the cylindrical electrodes.
        a (float): Radius of the inner electrode.
        b (float): Radius of the outer electrode.
    """
    if a >= b:
        print("Radius 'a' must be less than 'b'.")
        return

    # Common denominator term
    denominator_factor = 1 - (a**2 / b**2)
    
    # Check for division by zero
    if denominator_factor == 0:
        print("The denominator (1 - a^2/b^2) cannot be zero.")
        return

    # Total volume charge from Option B
    q_v = (-4 * V * epsilon * math.pi * L) / denominator_factor
    
    # Total surface charge on the inner electrode from Option B
    q_s_a = (2 * math.pi * L * V * epsilon) / denominator_factor

    # Total surface charge on the outer electrode from Option B
    q_s_b = (-4 * math.pi * L * V * epsilon * a**2) / (b**2 * denominator_factor)

    print("Based on Option B:")
    print(f"Total volume charge (q_v) = {q_v:.4f}")
    print(f"  Formula: q_v = -4 * V * epsilon * pi * L / (1 - a^2/b^2)")
    print(f"  = -4 * {V} * {epsilon} * pi * {L} / (1 - {a}^2/{b}^2)")
    
    print(f"Total surface charge on inner electrode (q_s(a)) = {q_s_a:.4f}")
    print(f"  Formula: q_s(r=a) = 2 * pi * L * V * epsilon / (1 - a^2/b^2)")
    print(f"  = 2 * pi * {L} * {V} * {epsilon} / (1 - {a}^2/{b}^2)")

    print(f"Total surface charge on outer electrode (q_s(b)) = {q_s_b:.4f}")
    print(f"  Formula: q_s(r=b) = -4 * pi * L * V * epsilon * a^2 / (b^2 * (1 - a^2/b^2))")
    print(f"  = -4 * pi * {L} * {V} * {epsilon} * {a}^2 / ({b}^2 * (1 - {a}^2/{b}^2))")

# Example values for demonstration
V_val = 100.0
epsilon_val = 8.854e-12  # Permittivity of free space as an example
L_val = 0.5
a_val = 0.1
b_val = 0.2

calculate_charges(V_val, epsilon_val, L_val, a_val, b_val)