import math

def calculate_charges(V, epsilon, L, a, b):
    """
    Calculates the total free volume charge and surface charges on electrodes
    based on the formulas in the selected answer choice B.

    Args:
        V (float): Applied DC voltage.
        epsilon (float): Permittivity of the Ohmic material.
        L (float): Depth of the cylindrical electrodes.
        a (float): Radius of the inner electrode.
        b (float): Radius of the outer electrode.
    """
    if not (b > a > 0 and L > 0):
        print("Invalid dimensions. Please ensure b > a > 0 and L > 0.")
        return

    # Common factors
    pi = math.pi
    factor_denom = 1 - (a**2 / b**2)

    # --- Calculations based on Option B ---

    # Total volume charge (q_v)
    # Equation: q_v = -4 * pi * L * V * epsilon / (1 - a^2/b^2)
    qv_numerator = -4 * pi * L * V * epsilon
    q_v = qv_numerator / factor_denom
    print("Total volume charge (q_v):")
    print(f"Equation: q_v = (-4 * pi * L * V * epsilon) / (1 - a^2/b^2)")
    print(f"q_v = ({qv_numerator:.4e}) / ({factor_denom:.4f}) = {q_v:.4e} C\n")

    # Total surface charge on the inner electrode (q_s at r=a)
    # Equation: q_s(a) = 2 * pi * L * V * epsilon / (1 - a^2/b^2)
    qsa_numerator = 2 * pi * L * V * epsilon
    q_s_a = qsa_numerator / factor_denom
    print("Total surface charge on inner electrode (q_s(r=a)):")
    print(f"Equation: q_s(a) = (2 * pi * L * V * epsilon) / (1 - a^2/b^2)")
    print(f"q_s(a) = ({qsa_numerator:.4e}) / ({factor_denom:.4f}) = {q_s_a:.4e} C\n")

    # Total surface charge on the outer electrode (q_s at r=b)
    # Equation: q_s(b) = -4 * pi * L * V * epsilon * a^2 / (b^2 * (1 - a^2/b^2))
    qsb_numerator = -4 * pi * L * V * epsilon * a**2
    qsb_denominator = b**2 * factor_denom
    q_s_b = qsb_numerator / qsb_denominator
    print("Total surface charge on outer electrode (q_s(r=b)):")
    print(f"Equation: q_s(b) = (-4 * pi * L * V * epsilon * a^2) / (b^2 * (1 - a^2/b^2))")
    print(f"q_s(b) = ({qsb_numerator:.4e}) / ({qsb_denominator:.4e}) = {q_s_b:.4e} C\n")


# --- Example Usage ---
# You can change these values to see the results for a specific case.
# Using SI units.
VOLTAGE = 100.0  # Volts
PERMITTIVITY = 2.5 * 8.854e-12  # F/m (relative permittivity of 2.5)
LENGTH = 0.5  # meters
RADIUS_A = 0.01  # meters (1 cm)
RADIUS_B = 0.03  # meters (3 cm)

print("Calculating charges with the following parameters:")
print(f"V = {VOLTAGE} V")
print(f"epsilon = {PERMITTIVITY:.4e} F/m")
print(f"L = {LENGTH} m")
print(f"a = {RADIUS_A} m")
print(f"b = {RADIUS_B} m\n")

calculate_charges(VOLTAGE, PERMITTIVITY, LENGTH, RADIUS_A, RADIUS_B)