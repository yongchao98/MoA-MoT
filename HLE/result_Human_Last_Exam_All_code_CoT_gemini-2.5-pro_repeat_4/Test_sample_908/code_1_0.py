import math

def calculate_charges(V, eps, L, a, b):
    """
    Calculates the volume and surface charges based on the formulas in Option B.

    Args:
        V (float): DC voltage in Volts.
        eps (float): Permittivity of the material in F/m.
        L (float): Length of the cylinders in meters.
        a (float): Radius of the inner cylinder in meters.
        b (float): Radius of the outer cylinder in meters.
    """
    pi = math.pi
    
    # Common denominator term
    denominator = 1 - (a**2 / b**2)
    
    if denominator == 0:
        print("Error: Radii a and b cannot be equal.")
        return

    # Calculate total volume charge based on option B
    q_v = (-4 * V * eps * pi * L) / denominator
    
    # Calculate total surface charge on inner electrode based on option B
    q_s_a = (2 * pi * L * V * eps) / denominator
    
    # Calculate total surface charge on outer electrode based on option B
    q_s_b = (-4 * pi * L * V * eps * a**2) / (b**2 * denominator)
    
    print("Assuming Option B is the intended answer:")
    print(f"Total volume charge (q_v) = (-4 * {V} * {eps:.3e} * {pi:.4f} * {L}) / (1 - {a}^2/{b}^2) = {q_v:.3e} C")
    print(f"Total surface charge on inner electrode (q_s(r=a)) = (2 * {pi:.4f} * {L} * {V} * {eps:.3e}) / (1 - {a}^2/{b}^2) = {q_s_a:.3e} C")
    print(f"Total surface charge on outer electrode (q_s(r=b)) = (-4 * {pi:.4f} * {L} * {V} * {eps:.3e} * {a}^2) / ({b}^2 * (1 - {a}^2/{b}^2)) = {q_s_b:.3e} C")

# Example usage with placeholder values
V_val = 10.0  # Volts
eps_r = 2.0  # Relative permittivity
eps_0 = 8.854e-12 # Permittivity of free space F/m
eps_val = eps_r * eps_0
L_val = 0.1  # meters
a_val = 0.01 # meters
b_val = 0.02 # meters

calculate_charges(V_val, eps_val, L_val, a_val, b_val)