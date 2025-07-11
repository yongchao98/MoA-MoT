import math

def calculate_reactive_power():
    """
    Calculates the reactive power needed for voltage compensation at Bus 4.
    """
    # Step 1: Define the given values from the problem description.
    V_kV = 400  # Target bus voltage in kilovolts (kV)
    Xc = 5      # Reactance of the shunt capacitor bank in Ohms

    # Step 2: Convert the voltage from kilovolts (kV) to volts (V) for the calculation.
    V_V = V_kV * 1000

    # Step 3: Calculate the reactive power (Qc) supplied by the capacitor bank using the formula Qc = V^2 / Xc.
    # This formula uses the line-to-line voltage.
    Qc_VAR = V_V**2 / Xc

    # Step 4: Convert the result from VAR to Megavars (MVAR) for a more standard representation.
    Qc_MVAR = Qc_VAR / 1e6

    # Step 5: Print the explanation, the final equation with values, and the result.
    print("The reactive power (Qc) supplied by the shunt capacitor is calculated to maintain the voltage at Bus 4.")
    print("The formula used is: Qc = V^2 / Xc\n")
    
    print("Given values:")
    print(f"  - Target Voltage (V) = {V_kV} kV")
    print(f"  - Capacitor Reactance (Xc) = {Xc} Ω\n")

    print("Calculation steps:")
    print(f"1. Convert voltage to Volts: V = {V_kV} kV * 1000 = {V_V:,.0f} V")
    print(f"2. Substitute values into the equation:")
    print(f"   Qc = ({V_V:,.0f} V)^2 / {Xc} Ω")
    print(f"   Qc = {V_V**2:,.0f} / {Xc}")
    print(f"   Qc = {Qc_VAR:,.0f} VAR\n")
    
    print("3. Convert VAR to MVAR:")
    print(f"   Qc = {Qc_VAR:,.0f} VAR / 1,000,000 = {Qc_MVAR:,.0f} MVAR\n")

    print(f"Therefore, the reactive power needed from the capacitor bank is {Qc_MVAR:,.0f} MVAR.")
    
    # Output the final answer in the specified format.
    print(f"\n<<<{Qc_MVAR:,.0f}>>>")

# Execute the function
calculate_reactive_power()