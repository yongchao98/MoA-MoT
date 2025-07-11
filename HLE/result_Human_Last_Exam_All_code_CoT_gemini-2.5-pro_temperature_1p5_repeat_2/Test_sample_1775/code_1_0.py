def calculate_displacement_field():
    """
    Calculates the electric displacement field in a dual-gate FET.
    
    This script prompts the user for the top and back gate capacitances per area (Ctg, Cbg)
    and the applied gate voltages (Vtg, Vbg) to calculate the net displacement field (D)
    through the grounded transistor channel.
    
    The formula used is D = C_bg * V_bg - C_tg * V_tg, where the positive direction
    is defined from the back gate towards the top gate.
    """
    try:
        # Get user input for the parameters
        V_tg_str = input("Enter the top gate voltage (V_tg) in Volts: ")
        V_tg = float(V_tg_str)
        
        C_tg_str = input("Enter the top gate capacitance per area (C_tg) in F/m^2: ")
        C_tg = float(C_tg_str)
        
        V_bg_str = input("Enter the back gate voltage (V_bg) in Volts: ")
        V_bg = float(V_bg_str)
        
        C_bg_str = input("Enter the back gate capacitance per area (C_bg) in F/m^2: ")
        C_bg = float(C_bg_str)
        
        # Calculate the displacement field
        # D = C_bg * V_bg - C_tg * V_tg
        displacement_field = C_bg * V_bg - C_tg * V_tg
        
        # Output the results
        print("\n--- Calculation ---")
        print("The formula for the displacement field (D) is:")
        print("D = (Back Gate Capacitance * Back Gate Voltage) - (Top Gate Capacitance * Top Gate Voltage)")
        print("D = C_bg * V_bg - C_tg * V_tg\n")
        
        print("Substituting the given values:")
        # Print the final equation with the numbers plugged in
        print(f"D = ({C_bg} F/m^2 * {V_bg} V) - ({C_tg} F/m^2 * {V_tg} V)")
        
        print("\n--- Result ---")
        print(f"The net displacement field (D) is: {displacement_field:.4e} C/m^2")
        print("(A positive value indicates the field points from the back gate to the top gate)")

    except ValueError:
        print("\nError: Invalid input. Please enter valid numbers.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == "__main__":
    calculate_displacement_field()