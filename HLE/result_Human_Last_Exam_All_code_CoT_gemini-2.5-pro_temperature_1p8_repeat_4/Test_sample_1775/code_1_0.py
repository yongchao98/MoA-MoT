import sys

def calculate_displacement_field(c_tg, c_bg, v_tg, v_bg):
    """
    Calculates the displacement field in a dual-gate FET using a thin-channel model.

    Args:
        c_tg (float): Top gate capacitance per unit area (in F/m^2).
        c_bg (float): Back gate capacitance per unit area (in F/m^2).
        v_tg (float): Top gate voltage (in V).
        v_bg (float): Back gate voltage (in V).

    Returns:
        float: The displacement field D (in C/m^2).
    """
    # The displacement field D is given by the formula:
    # D = (C_tg * C_bg) / (C_tg + C_bg) * (V_tg - V_bg)
    # This formula can be derived by modeling the transistor as a thin, neutral floating channel
    # whose potential is determined by capacitive coupling to the gates.

    if (c_tg + c_bg) == 0:
        print("Error: Sum of capacitances cannot be zero.", file=sys.stderr)
        return None

    # Calculate the combined capacitance term
    c_eff = (c_tg * c_bg) / (c_tg + c_bg)
    
    # Calculate the voltage difference
    delta_v = v_tg - v_bg
    
    # Calculate the displacement field
    d_field = c_eff * delta_v
    
    return d_field

# --- Example Usage ---
# Define example parameters for the transistor.
# Capacitances are often given in uF/cm^2, we convert them to F/m^2 for SI units.
# 1 uF/cm^2 = 1e-6 F / (1e-2 m)^2 = 1e-6 / 1e-4 F/m^2 = 1e-2 F/m^2.
c_tg_val = 1.5e-2  # Corresponds to 1.5 uF/cm^2
c_bg_val = 0.5e-2  # Corresponds to 0.5 uF/cm^2
v_tg_val = 1.0     # 1.0 Volts
v_bg_val = -2.0    # -2.0 Volts

# Calculate the displacement field
displacement_field = calculate_displacement_field(c_tg_val, c_bg_val, v_tg_val, v_bg_val)

if displacement_field is not None:
    # Print the final result in a descriptive equation
    print("The displacement field (D) is calculated based on the thin-channel approximation.")
    print("Formula: D = (C_tg * C_bg) / (C_tg + C_bg) * (V_tg - V_bg)\n")
    print("Calculation with provided values:")
    # Using f-string formatting to display the equation with the numbers plugged in
    print(f"D = ({c_tg_val:.3e} F/m^2 * {c_bg_val:.3e} F/m^2) / ({c_tg_val:.3e} F/m^2 + {c_bg_val:.3e} F/m^2) * ({v_tg_val:.1f} V - ({v_bg_val:.1f} V))")
    
    c_sum = c_tg_val + c_bg_val
    c_prod = c_tg_val * c_bg_val
    delta_v = v_tg_val - v_bg_val
    
    print(f"D = ({c_prod:.3e}) / ({c_sum:.3e}) * ({delta_v:.1f})")
    print(f"D = {displacement_field:.4e} C/m^2")

<<<1.1250e-02>>>