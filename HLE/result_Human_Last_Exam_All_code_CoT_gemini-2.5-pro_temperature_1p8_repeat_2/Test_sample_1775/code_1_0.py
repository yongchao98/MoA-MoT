import sys

def calculate_displacement_field():
    """
    Calculates the net displacement field in a dual-gate FET.

    The function models the transistor as a grounded channel between a top gate and
    a back gate. The net displacement field is the superposition of the fields
    from the top and back gates.

    We define the direction from the back gate to the top gate as positive.
    - The field from the top gate (V_tg > 0) points down (negative direction).
    - The field from the back gate (V_bg > 0) points up (positive direction).

    The formula used is: D_net = (C_bg * V_bg) - (C_tg * V_tg).

    The dielectric constant of the transistor (epsilon_s) is not required for
    this model, as the net external displacement field is determined by the
    gate voltages and capacitances per area.
    """
    try:
        # Prompting for user input with default values for demonstration
        print("Enter the values for the transistor parameters (leave blank for default values).")
        
        ctg_str = input("Enter top gate capacitance per area (C_tg) in F/m^2 [default: 2e-6]: ")
        C_tg = float(ctg_str) if ctg_str else 2e-6

        vtg_str = input("Enter top gate voltage (V_tg) in V [default: 3.0]: ")
        V_tg = float(vtg_str) if vtg_str else 3.0

        cbg_str = input("Enter back gate capacitance per area (C_bg) in F/m^2 [default: 1.5e-6]: ")
        C_bg = float(cbg_str) if cbg_str else 1.5e-6

        vbg_str = input("Enter back gate voltage (V_bg) in V [default: 2.0]: ")
        V_bg = float(vbg_str) if vbg_str else 2.0

        # Calculate the individual components
        D_top_component = C_tg * V_tg
        D_bot_component = C_bg * V_bg

        # Calculate the net displacement field
        D_net = D_bot_component - D_top_component

        # --- Output the results and the final equation ---
        print("\n--- Calculation ---")
        print("The net displacement field (D_net) is calculated as the sum of the fields from each gate.")
        print("We define the direction from the back gate to the top gate as positive.")
        print("\nFormula: D_net = (C_bg * V_bg) - (C_tg * V_tg)")
        
        # Print the final equation with all numbers
        print("\nStep-by-step substitution:")
        print(f"D_net = ({C_bg} * {V_bg}) - ({C_tg} * {V_tg})")
        print(f"D_net = {D_bot_component} - {D_top_component}")
        print(f"D_net = {D_net:.2e} C/m^2") # Units of displacement field are Coulombs per square meter

        print("\n--- Conclusion ---")
        if D_net > 0:
            print(f"The final displacement field is {D_net:.2e} C/m^2.")
            print("The positive sign indicates the net field points from the back gate towards the top gate.")
        elif D_net < 0:
            print(f"The final displacement field is {D_net:.2e} C/m^2.")
            print("The negative sign indicates the net field points from the top gate towards the back gate.")
        else:
            print("The final displacement field is 0.0 C/m^2.")
            print("The fields from the top and back gates are perfectly balanced.")

    except ValueError:
        print("\nError: Invalid input. Please enter numeric values.", file=sys.stderr)
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}", file=sys.stderr)

if __name__ == '__main__':
    calculate_displacement_field()
