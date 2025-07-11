import math

def calculate_power_loss():
    """
    Calculates the total resistive power losses in the transmission lines
    based on the power system's specifications.
    """
    # --- Step 1: Define Constants from the problem ---

    # Generator parameters
    N_gen = 3
    S_gen_rated_per_unit = 100.0  # MVA
    PF = 0.9  # lagging

    # Load parameters (Real Power)
    P_base_load_per_unit = 50.0 / 3.0  # MW
    N_base_loads = 3
    P_add_load = 100.0  # MW

    # Transmission line resistances from the table in Ohms
    R_lines = {
        'T7_8': 4.50,
        'T8_9': 6.30,
        'T5_7': 16.93,
        'T6_9': 20.63,
        'T4_5': 5.29,
        'T4_6': 8.99,
    }

    # Transformer parameters
    N_tx = 3
    S_tx_base = 210.0  # MVA
    V_tx_base_hv = 230.0  # kV
    Z_tx_percent = 8.0  # %
    # Assumption for transformer X/R ratio, as it's not provided.
    # A value of 15 is typical for large power transformers.
    X_div_R_ratio_tx = 15.0

    # --- Step 2: Calculate Total Generated Power ---
    P_gen_per_unit = S_gen_rated_per_unit * PF
    P_gen_total = N_gen * P_gen_per_unit

    # --- Step 3: Calculate Total Load Power ---
    P_load_total = (N_base_loads * P_base_load_per_unit) + P_add_load

    # --- Step 4: Calculate Total System Power Loss ---
    P_loss_total = P_gen_total - P_load_total

    # --- Step 5 & 6: Calculate Total Line Resistance ---
    R_lines_total = sum(R_lines.values())

    # --- Step 7: Calculate Single Transformer Resistance ---
    # Base impedance Z_base = V_base^2 / S_base
    # Note: (kV^2 / MVA) directly gives Ohms
    Z_base_tx = (V_tx_base_hv ** 2) / S_tx_base
    # Actual impedance Z_actual = (%Z/100) * Z_base
    Z_actual_tx = (Z_tx_percent / 100.0) * Z_base_tx
    # From Z^2 = R^2 + X^2 and X = k*R (where k=X/R), we get R = sqrt(Z^2 / (1 + k^2))
    R_tx_one = math.sqrt(Z_actual_tx**2 / (1 + X_div_R_ratio_tx**2))

    # --- Step 8: Calculate Total Transformer Resistance ---
    R_tx_total = N_tx * R_tx_one

    # --- Step 9: Apportion Loss to find Line Loss ---
    # P_loss_lines = P_loss_total * (R_lines_total / (R_lines_total + R_tx_total))
    P_loss_lines = P_loss_total * (R_lines_total / (R_lines_total + R_tx_total))
    
    # --- Step 10: Print the final calculation and result ---
    print("Calculation of Total Resistive Power Losses in Transmission Lines")
    print("-" * 60)
    print(f"Total Generated Power (P_gen_total): {P_gen_total:.3f} MW")
    print(f"Total Load Power (P_load_total): {P_load_total:.3f} MW")
    print(f"Total System Real Power Loss (P_loss_total) = {P_gen_total:.3f} - {P_load_total:.3f} = {P_loss_total:.3f} MW")
    print("-" * 60)
    print("Component Resistance Estimation:")
    print(f"Total Line Resistance (R_lines_total): {R_lines_total:.3f} Ω")
    print(f"Total Transformer Resistance (R_tx_total): {R_tx_total:.3f} Ω (assuming an X/R ratio of {X_div_R_ratio_tx})")
    print(f"Total System Resistance (R_lines_total + R_tx_total): {R_lines_total + R_tx_total:.3f} Ω")
    print("-" * 60)
    
    # --- Final Equation Output ---
    print("The final calculation is performed as follows:")
    print("P_loss_lines = P_loss_total * (R_lines_total / (R_lines_total + R_tx_total))")
    print(f"P_loss_lines = {P_loss_total:.3f} * ({R_lines_total:.3f} / ({R_lines_total:.3f} + {R_tx_total:.3f}))")
    final_ratio = R_lines_total / (R_lines_total + R_tx_total)
    print(f"P_loss_lines = {P_loss_total:.3f} * {final_ratio:.3f}")

    print("\nResult:")
    print(f"The total power loss in the transmission lines is: {P_loss_lines:.3f} MW")

if __name__ == '__main__':
    calculate_power_loss()