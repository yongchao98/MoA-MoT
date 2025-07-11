import math

def calculate_charges(V, epsilon, L, a, b):
    """
    Calculates the total volume and surface charges based on the formulas from option B.

    Args:
        V (float): DC voltage applied.
        epsilon (float): Permittivity of the Ohmic material.
        L (float): Depth of the electrodes.
        a (float): Radius of the inner electrode.
        b (float): Radius of the outer electrode.
    """
    pi = math.pi
    
    # Check for valid radii
    if a >= b:
        print("Error: Radius 'a' must be smaller than radius 'b'.")
        return

    # Denominator term
    denominator = 1 - (a**2 / b**2)
    
    # --- Calculate Total Volume Charge (q_v) from Option B ---
    q_v_numerator = -4 * V * epsilon * pi * L
    q_v = q_v_numerator / denominator
    print("--- Option B Calculations ---")
    print(f"Total volume charge = q_v = (-4 * {V} * {epsilon:.3e} * {pi:.4f} * {L}) / (1 - ({a}^2 / {b}^2))")
    print(f"q_v = {q_v_numerator:.4e} / {denominator:.4f} = {q_v:.4e} C\n")
    
    # --- Calculate Total Surface Charge on Inner Electrode (q_s(a)) from Option B ---
    qs_a_numerator = 2 * pi * L * V * epsilon
    qs_a = qs_a_numerator / denominator
    print(f"Total surface charge on inner electrode = q_s(a) = (2 * {pi:.4f} * {L} * {V} * {epsilon:.3e}) / (1 - ({a}^2 / {b}^2))")
    print(f"q_s(a) = {qs_a_numerator:.4e} / {denominator:.4f} = {qs_a:.4e} C\n")

    # --- Calculate Total Surface Charge on Outer Electrode (q_s(b)) from Option B ---
    qs_b_numerator = -4 * pi * L * V * epsilon * (a**2)
    qs_b_denominator = (b**2) * denominator
    qs_b = qs_b_numerator / qs_b_denominator
    print(f"Total surface charge on outer electrode = q_s(b) = (-4 * {pi:.4f} * {L} * {V} * {epsilon:.3e} * {a}^2) / ({b}^2 * (1 - ({a}^2 / {b}^2)))")
    print(f"q_s(b) = {qs_b_numerator:.4e} / {qs_b_denominator:.4f} = {qs_b:.4e} C\n")


# --- Example Usage ---
# Define placeholder values for the variables
V_val = 100.0  # Volts
epsilon_val = 2.5 * 8.854e-12  # F/m (relative permittivity of 2.5)
L_val = 0.5  # meters
a_val = 0.01 # meters
b_val = 0.02 # meters

# Run the calculation and print the results
calculate_charges(V_val, epsilon_val, L_val, a_val, b_val)