def calculate_displacement_field():
    """
    Calculates the total displacement field in a dual-gate FET.
    """
    print("This script calculates the net displacement field (D) in a dual-gate transistor.")
    print("The system is modeled as two capacitors connected to the grounded channel.")
    print("The formula is: D_total = (C_tg * V_tg) + (C_bg * V_bg)\n")

    try:
        # Prompt user for input values
        c_tg_str = input("Enter the top gate capacitance per area (C_tg in F/m^2): ")
        v_tg_str = input("Enter the top gate voltage (V_tg in Volts): ")
        c_bg_str = input("Enter the back gate capacitance per area (C_bg in F/m^2): ")
        v_bg_str = input("Enter the back gate voltage (V_bg in Volts): ")

        # Convert input strings to floating-point numbers
        c_tg = float(c_tg_str)
        v_tg = float(v_tg_str)
        c_bg = float(c_bg_str)
        v_bg = float(v_bg_str)

        # Perform the calculation
        d_tg = c_tg * v_tg
        d_bg = c_bg * v_bg
        d_total = d_tg + d_bg

        # Display the equation with the user's numbers
        print("\n--- Calculation Breakdown ---")
        print(f"The final equation is given by D = (C_tg * V_tg) + (C_bg * V_bg)")
        print(f"Plugging in the values: D = ({c_tg} * {v_tg}) + ({c_bg} * {v_bg})")
        print(f"Contribution from each gate: D = {d_tg} + {d_bg}")

        # Display the final result
        print("\n--- Final Result ---")
        print(f"The total displacement field (D) is: {d_total} C/m^2")

    except ValueError:
        print("\nError: Invalid input. Please make sure to enter numerical values only.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == "__main__":
    calculate_displacement_field()