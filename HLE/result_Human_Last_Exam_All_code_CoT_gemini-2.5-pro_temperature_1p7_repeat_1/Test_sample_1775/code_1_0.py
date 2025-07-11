import sys

def calculate_displacement_field():
    """
    Calculates the displacement field in a dual-gate FET.
    """
    try:
        # Prompt user for input values
        print("Please provide the following values for the transistor:")
        Ctg = float(input("Enter the top gate capacitance per area (C_tg) in F/m^2: "))
        Vtg = float(input("Enter the top gate voltage (V_tg) in V: "))
        Cbg = float(input("Enter the back gate capacitance per area (C_bg) in F/m^2: "))
        Vbg = float(input("Enter the back gate voltage (V_bg) in V: "))

        # --- Calculation ---
        # Displacement field from the top gate
        D_tg = Ctg * Vtg
        # Displacement field from the back gate
        D_bg = Cbg * Vbg
        # Total displacement field by superposition
        D_total = D_tg + D_bg

        # --- Output the results ---
        print("\n-------------------------------------------------")
        print("Calculating the Total Displacement Field (D)")
        print("-------------------------------------------------")
        print("The displacement field is the sum of the contributions from the top and back gates.")
        print("Formula: D_total = C_tg * V_tg + C_bg * V_bg")
        
        print("\nSubstituting the given values into the equation:")
        # Print the equation with the user's numbers
        print(f"D_total = ({Ctg} F/m^2) * ({Vtg} V) + ({Cbg} F/m^2) * ({Vbg} V)")
        
        print("\nCalculating each term:")
        print(f"Contribution from top gate: {D_tg:.4g} C/m^2")
        print(f"Contribution from back gate: {D_bg:.4g} C/m^2")
        
        # Print the final combined equation and result
        print(f"\nD_total = {D_tg:.4g} C/m^2 + {D_bg:.4g} C/m^2")
        print(f"\nThe total displacement field is: {D_total:.4g} C/m^2")
        print("-------------------------------------------------\n")

    except ValueError:
        print("\nError: Invalid input. Please enter numerical values.", file=sys.stderr)
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}", file=sys.stderr)

if __name__ == '__main__':
    calculate_displacement_field()
