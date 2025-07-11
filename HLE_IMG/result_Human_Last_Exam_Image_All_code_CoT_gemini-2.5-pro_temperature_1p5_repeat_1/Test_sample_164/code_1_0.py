import math

def calculate_reactive_power():
    """
    Calculates the reactive power supplied by a shunt capacitor bank.
    
    The problem asks for the reactive power needed at Bus 4 to maintain its voltage
    at 400 kV, with the compensation provided by a shunt capacitor bank with a
    reactance of 5 Ohms. The most direct calculation is to determine the reactive
    power this specific capacitor will supply at the target voltage.
    """
    
    # --- Given Parameters ---
    V_kv = 400.0  # Target line-to-line voltage in kilovolts (kV)
    Xc_ohm = 5.0  # Reactance of the shunt capacitor bank in Ohms

    # --- Calculation ---
    # The formula for total three-phase reactive power in MVAR can be expressed as:
    # Qc (MVAR) = (V_ll_kV)^2 / Xc_ohm
    # This works because V_ll_kV^2 gives units of (kV)^2, and we need to convert
    # to MVAR. Qc = (V_ll)^2 / Xc = (V_kV * 1000)^2 / Xc = V_kV^2 * 1e6 / Xc.
    # To get MVAR, we divide by 1e6, so Qc_MVAR = V_kV^2 / Xc.
    
    Qc_mvar = V_kv**2 / Xc_ohm
    
    # --- Output Results ---
    print("Performing reactive power compensation analysis at Bus 4.")
    print(f"Target bus voltage: {int(V_kv)} kV")
    print(f"Shunt capacitor reactance (Xc): {int(Xc_ohm)} Ω")
    print("\nThe reactive power (Qc) supplied by the capacitor bank is calculated using the formula:")
    print("Qc (MVAR) = (Voltage in kV)^2 / (Reactance in Ω)")
    
    print("\nPlugging in the values:")
    # The final instruction asks to output each number in the final equation.
    # We will show the calculation step by step.
    V_squared = V_kv**2
    print(f"1. Square the voltage: {int(V_kv)}^2 = {int(V_squared)}")
    print(f"2. Divide by the reactance: {int(V_squared)} / {int(Xc_ohm)} = {int(Qc_mvar)}")
    
    print("\nThe final equation with numbers is:")
    print(f"{int(Qc_mvar)} = {int(V_kv)}^2 / {int(Xc_ohm)}")

    print(f"\nTherefore, the reactive power supplied by the capacitor bank is {int(Qc_mvar)} MVAR.")

if __name__ == "__main__":
    calculate_reactive_power()
    print("\n<<<32000>>>")
