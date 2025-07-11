import sys
import io

# Define a function to encapsulate the logic
def calculate_displacement_field():
    # --- Input Parameters ---
    # Top gate capacitance per unit area in Farads per square meter (F/m^2)
    Ctg = 1.725e-2

    # Top gate voltage in Volts (V)
    Vtg = 2.5

    # Back gate capacitance per unit area in Farads per square meter (F/m^2)
    Cbg = 0.345e-2

    # Back gate voltage in Volts (V)
    Vbg = 10.0
    # --------------------------

    # Step 1: Calculate the displacement field component from the top gate
    D_tg = Ctg * Vtg

    # Step 2: Calculate the displacement field component from the back gate
    D_bg = Cbg * Vbg

    # Step 3: Sum the components to get the total displacement field
    D_total = D_tg + D_bg

    # --- Output ---
    # Use a string buffer to capture the output
    buffer = io.StringIO()
    # Redirect stdout to the buffer
    original_stdout = sys.stdout
    sys.stdout = buffer

    print("The displacement field (D) through the transistor is the sum of the fields induced by the top and back gates.")
    print("The formula is: D = C_tg * V_tg + C_bg * V_bg")
    print("\nPlugging in the values:")
    print(f"D = ({Ctg}) * ({Vtg}) + ({Cbg}) * ({Vbg})")
    print(f"D = {D_tg} + {D_bg}")
    print(f"\nThe total displacement field is:")
    print(f"D = {D_total} C/m^2")
    
    # Restore stdout
    sys.stdout = original_stdout
    # Get the string from the buffer
    output_string = buffer.getvalue()
    
    # Print the string to the actual console
    print(output_string)

    # Return the final numerical answer in the required format
    return D_total

# Execute the function and capture the final answer
final_answer = calculate_displacement_field()

# Output the final answer in the specified format
# The value is returned as a string to handle scientific notation if necessary.
print(f"<<<{final_answer}>>>")
