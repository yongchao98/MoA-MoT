import math

def calculate_charges(V, eps, L, a, b):
    """
    Calculates the charges based on the formulas in option B.
    Note: These formulas are inconsistent with a first-principles physics derivation,
    which suggests errors in the problem's options.
    """
    pi = math.pi
    
    # Common denominator term
    denom = 1 - (a**2 / b**2)
    
    # Formulas from Option B
    q_v = (-4 * V * eps * pi * L) / denom
    q_s_a = (2 * pi * L * V * eps) / denom
    q_s_b = (-4 * pi * L * V * eps * a**2) / (b**2 * denom)

    print("Calculations are based on the formulas provided in Option B.")
    print("Total volume charge q_v = (-4 * pi * L * V * eps) / (1 - (a/b)**2)")
    print(f"q_v = {q_v:.4e} C\n")
    
    print("Total surface charge on inner electrode q_s(a) = (2 * pi * L * V * eps) / (1 - (a/b)**2)")
    print(f"q_s(r=a) = {q_s_a:.4e} C\n")
    
    print("Total surface charge on outer electrode q_s(b) = (-4 * pi * L * V * eps * a**2) / (b**2 * (1 - (a/b)**2))")
    print(f"q_s(r=b) = {q_s_b:.4e} C\n")

# Example parameters for calculation
V_val = 100.0  # Volts
eps_0 = 8.854e-12 # F/m
eps_r = 2.0 # Relative permittivity
eps_val = eps_r * eps_0
L_val = 0.5  # meters
a_val = 0.01 # meters
b_val = 0.02 # meters

calculate_charges(V_val, eps_val, L_val, a_val, b_val)