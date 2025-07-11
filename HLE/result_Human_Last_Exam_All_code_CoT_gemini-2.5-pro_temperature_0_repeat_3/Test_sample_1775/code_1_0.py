import numpy as np

def calculate_displacement_field():
    """
    Calculates the net displacement field in a dual-gate FET.

    The function uses a physical model where the net displacement field (D)
    is the superposition of the fields created by the top and back gates.
    The formula used is D = C_bg * V_bg - C_tg * V_tg, where the positive
    direction is defined from the back gate towards the top gate.
    """

    # --- User-defined variables (example values) ---
    # Top gate voltage in Volts (V)
    V_tg = 1.2
    # Back gate voltage in Volts (V)
    V_bg = 0.8
    # Top gate capacitance per area in Farads per square meter (F/m^2)
    # Example: 3.9 * 8.854e-12 F/m / 5e-9 m = 6.9e-3 F/m^2
    C_tg = 6.9e-3
    # Back gate capacitance per area in Farads per square meter (F/m^2)
    # Example: 3.9 * 8.854e-12 F/m / 10e-9 m = 3.45e-3 F/m^2
    C_bg = 3.45e-3

    # The dielectric constant of the transistor material (epsilon_s) is not
    # needed for this calculation, as the displacement field D is continuous
    # across the dielectric-semiconductor boundary.

    # --- Calculation ---
    # Calculate the net displacement field using the superposition principle.
    D_net = C_bg * V_bg - C_tg * V_tg

    # --- Output ---
    print("To find the displacement field (D) through the transistor, we use the principle of superposition.")
    print("The net field is the sum of the fields from the back gate and the top gate.")
    print("Defining the direction from back gate to top gate as positive, the formula is:")
    print("D = (C_bg * V_bg) - (C_tg * V_tg)\n")

    print("Given the following values:")
    print(f"  Top Gate Voltage (V_tg):    {V_tg} V")
    print(f"  Back Gate Voltage (V_bg):   {V_bg} V")
    print(f"  Top Gate Capacitance (C_tg):  {C_tg:.2e} F/m^2")
    print(f"  Back Gate Capacitance (C_bg): {C_bg:.2e} F/m^2\n")

    print("Plugging the numbers into the equation:")
    # This print statement shows the final equation with all the numbers, as requested.
    print(f"D = ({C_bg:.2e} F/m^2 * {V_bg} V) - ({C_tg:.2e} F/m^2 * {V_tg} V)")

    # Calculate the individual components for clarity
    D_bg_val = C_bg * V_bg
    D_tg_val = C_tg * V_tg
    print(f"D = ({D_bg_val:.3e}) - ({D_tg_val:.3e}) C/m^2")

    print(f"\nThe final displacement field is:")
    print(f"D = {D_net:.3e} C/m^2")

    print("\nNote: A positive result indicates the net field vector points from the back gate to the top gate.")
    print("A negative result indicates the net field vector points from the top gate to the back gate.")

if __name__ == '__main__':
    calculate_displacement_field()