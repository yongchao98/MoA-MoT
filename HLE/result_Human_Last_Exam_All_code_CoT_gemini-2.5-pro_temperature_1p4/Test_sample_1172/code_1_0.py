import numpy as np

def calculate_inductance_change(h, d, R1, R2):
    """
    Calculates the change in mutual inductance per unit length between two circuits
    when placed inside a magnetic concentrator.

    Args:
        h (float): Separation between wires in a circuit (meters).
        d (float): Distance between the two circuits (meters).
        R1 (float): Inner radius of the concentrator shell (meters).
        R2 (float): Outer radius of the concentrator shell (meters).
    """
    # Check for valid inputs
    if not (d > 0 and R1 > 0):
        print("Error: d and R1 must be positive.")
        return
    if d <= h:
        print("Warning: The calculation assumes d >> h.")
    if d >= R1:
        print("Warning: The calculation assumes both circuits are well inside the concentrator (d << R1).")

    # Physical constant: permeability of free space (H/m)
    mu_0 = 4 * np.pi * 1e-7

    # Step 1: Calculate the mutual inductance per unit length for the bare circuits (m1)
    # m1 = (μ₀ / 2π) * ln(1 + (h/d)²)
    m1_factor = mu_0 / (2 * np.pi)
    geom_factor = np.log(1 + (h / d)**2)
    m1 = m1_factor * geom_factor

    # Step 2: Calculate the change in mutual inductance (Δm)
    # Δm = ((R₂ - R₁) / R₁) * m1
    concentrator_factor = (R2 - R1) / R1
    delta_m = concentrator_factor * m1

    # Output the results, showing the numbers in the final equation
    print("--- Calculation of Mutual Inductance Change per Unit Length ---")
    print(f"Mutual inductance without concentrator (M1): {m1:.4e} H/m")
    print(f"Concentrator amplification factor ((R2 - R1) / R1): {concentrator_factor:.4f}")
    print(f"The change in mutual inductance (ΔM = M2 - M1) is:")
    print(f"ΔM = {concentrator_factor:.4f} * {m1:.4e} H/m")
    print(f"ΔM = {delta_m:.4e} H/m")
    
    # Returning the final value as a formatted string for the answer
    return f"((R2 - R1) / R1) * (mu_0 / (2*pi)) * ln(1 + (h/d)**2)"

# --- Example Usage ---
# You can change these values to match a specific scenario.
# All values are in meters.
h_val = 0.01   # 1 cm separation between wires
d_val = 0.1    # 10 cm separation between circuits
R1_val = 0.2   # 20 cm inner radius of concentrator
R2_val = 0.25  # 25 cm outer radius of concentrator

# Calculate and print the final expression and a numerical example
final_expression = calculate_inductance_change(h_val, d_val, R1_val, R2_val)
print("\nFinal symbolic expression for the change in mutual inductance per unit length:")
print(final_expression)

# The final answer format
final_answer_expression = "((R2 - R1) / R1) * (mu_0 / (2*pi)) * log(1 + (h/d)**2)"
# Let's adjust for the case where 'log' might be interpreted as log10
final_answer_expression_ln = "((R2 - R1) / R1) * (mu_0 / (2*pi)) * ln(1 + (h/d)**2)"

# The prompt seems to ask for the expression itself.
#<<<((R₂ - R₁) / R₁) * (μ₀ / (2π)) * ln(1 + (h/d)²))>>>
# My python code already produced this as a string, let's use that one.
# final_expression variable has: "((R2 - R1) / R1) * (mu_0 / (2*pi)) * ln(1 + (h/d)**2)"

# We have to be careful with the output format required by the problem
# "directly return the answer with the format <<<answer content>>>"
# The answer is the expression itself.
# Final Answer format should contain the final expression.

final_answer = f"<<<{final_answer_expression_ln}>>>"
# This will be appended at the end of the entire response.
