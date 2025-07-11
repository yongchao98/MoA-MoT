def calculate_displacement_field():
    """
    Calculates the total displacement field in a dual-gate FET.

    The displacement field D is calculated using the formula:
    D = C_tg * V_tg + C_bg * V_bg
    where:
    C_tg = Top gate capacitance per unit area (e.g., in F/m^2)
    V_tg = Top gate voltage (in V)
    C_bg = Back gate capacitance per unit area (e.g., in F/m^2)
    V_bg = Back gate voltage (in V)
    The transistor channel is assumed to be at ground potential (0V).
    """
    try:
        # Get input from the user
        ctg_str = input("Enter the top gate capacitance per area (C_tg): ")
        vtg_str = input("Enter the top gate voltage (V_tg): ")
        cbg_str = input("Enter the back gate capacitance per area (C_bg): ")
        vbg_str = input("Enter the back gate voltage (V_bg): ")

        # Convert inputs to floating-point numbers
        c_tg = float(ctg_str)
        v_tg = float(vtg_str)
        c_bg = float(cbg_str)
        v_bg = float(vbg_str)

        # Calculate the displacement field
        displacement_field = c_tg * v_tg + c_bg * v_bg

        # Print the equation with the plugged-in numbers and the final result
        # This shows each number in the final equation.
        print("\nThe displacement field D is calculated as:")
        print(f"D = (C_tg * V_tg) + (C_bg * V_bg)")
        print(f"D = ({c_tg} * {v_tg}) + ({c_bg} * {v_bg})")
        print(f"D = {c_tg * v_tg} + {c_bg * v_bg}")
        print(f"D = {displacement_field}")
        
        # The final answer in the specified format is the numerical result.
        # However, as per instructions, it's included as part of the thinking process
        # and the code calculates it. For the final output, just running the code suffices.
        
    except ValueError:
        print("\nError: Please enter valid numerical values.")
    except Exception as e:
        print(f"\nAn error occurred: {e}")

if __name__ == "__main__":
    calculate_displacement_field()