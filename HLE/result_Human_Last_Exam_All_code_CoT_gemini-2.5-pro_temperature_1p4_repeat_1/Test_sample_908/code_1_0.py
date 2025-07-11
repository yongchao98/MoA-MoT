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
    pi = math.pi
    
    # Common denominator term
    denominator_factor = 1 - (a**2 / b**2)
    
    if denominator_factor == 0:
        print("Error: a cannot be equal to b.")
        return

    # Total volume charge from Option B
    q_v = (-4 * V * epsilon * pi * L) / denominator_factor

    # Total surface charge on the inner electrode (r=a) from Option B
    q_s_a = (2 * pi * L * V * epsilon) / denominator_factor

    # Total surface charge on the outer electrode (r=b) from Option B
    q_s_b = (-4 * pi * L * V * epsilon * a**2) / (b**2 * denominator_factor)
    
    print("Based on the formulas in Answer Choice B:")
    print(f"Total volume charge (q_v) = {q_v}")
    print(f"Total surface charge on inner electrode (q_s(a)) = {q_s_a}")
    print(f"Total surface charge on outer electrode (q_s(b)) = {q_s_b}")

# Example usage with some placeholder values
V_val = 100.0  # Volts
epsilon_val = 8.854e-12 * 2.5 # Farads per meter (e.g., relative permittivity of 2.5)
L_val = 0.5   # meters
a_val = 0.01  # meters
b_val = 0.02  # meters

calculate_charges(V_val, epsilon_val, L_val, a_val, b_val)