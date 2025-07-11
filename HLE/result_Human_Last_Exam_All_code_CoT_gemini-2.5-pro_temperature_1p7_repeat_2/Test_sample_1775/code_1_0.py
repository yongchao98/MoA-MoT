import sys

def calculate_displacement_field(Ctg, Vtg, Cbg, Vbg):
    """
    Calculates the displacement field in a dual-gate FET.

    Args:
        Ctg (float): Top gate capacitance per area (e.g., in F/m^2).
        Vtg (float): Top gate voltage (in V).
        Cbg (float): Bottom gate capacitance per area (e.g., in F/m^2).
        Vbg (float): Bottom gate voltage (in V).
    """
    # The displacement field from the top gate
    D_tg = Ctg * Vtg
    # The displacement field from the bottom gate
    D_bg = Cbg * Vbg
    # The total displacement field is the sum of the fields terminating on the channel
    total_D = D_tg + D_bg

    # We print the result in the format: D = Ctg * Vtg + Cbg * Vbg = result
    # This shows the contribution of each term clearly.
    print(f"Displacement Field (D) = D_tg + D_bg")
    print(f"D = (Ctg * Vtg) + (Cbg * Vbg)")
    # Using 'g' format specifier to handle floating point and scientific notation nicely
    print(f"D = ({Ctg:g} * {Vtg:g}) + ({Cbg:g} * {Vbg:g}) = {total_D:g}")
    
    # The final answer to be parsed.
    # Note: Displacement field has units of charge per area (e.g., C/m^2)
    print(f"\nFinal Answer (in units of C/m^2 if inputs are in SI):")
    # Python cannot output a final format like <<<content>>>.
    # We will print the numerical value directly.
    # print(f'<<<{total_D:g}>>>')


# Example values based on common device parameters
# Capacitance per area is often given in uF/cm^2. 1 uF/cm^2 = 1e-2 F/m^2
# Let's use SI units (F/m^2) for calculation consistency.
Ctg_val = 2.5e-2 # F/m^2 (equivalent to 2.5 uF/cm^2)
Vtg_val = 1.5    # V
Cbg_val = 0.5e-2 # F/m^2 (equivalent to 0.5 uF/cm^2)
Vbg_val = 5.0    # V

# Execute the calculation with example values
# You can change these values to match a specific problem
calculate_displacement_field(Ctg_val, Vtg_val, Cbg_val, Vbg_val)

# The result of (2.5e-2 * 1.5) + (0.5e-2 * 5.0) is 0.0375 + 0.025 = 0.0625
final_answer = Ctg_val * Vtg_val + Cbg_val * Vbg_val
sys.stdout.write(f'<<<{final_answer:g}>>>')