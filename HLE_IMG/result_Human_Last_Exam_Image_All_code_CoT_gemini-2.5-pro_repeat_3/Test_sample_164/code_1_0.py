import math

def calculate_reactive_power_compensation():
    """
    Calculates the reactive power supplied by a shunt capacitor bank at Bus 4.
    
    The problem asks to determine the reactive power needed to maintain the voltage
    at Bus 4 at 400 kV, using a specified shunt capacitor bank. A full load flow
    analysis is not possible due to missing data (line lengths, load values).
    Therefore, the most direct interpretation is to calculate the reactive power
    supplied by the given capacitor at the target voltage.
    """
    
    # --- Given Parameters ---
    # Target line-to-line voltage at Bus 4 in Volts
    V_target_kv = 400.0
    V_target_v = V_target_kv * 1000
    
    # Reactance of the shunt capacitor bank in Ohms
    Xc = 5.0

    # --- Calculation ---
    # The formula for total 3-phase reactive power supplied by a shunt capacitor is:
    # Qc = V_ll^2 / Xc
    # where V_ll is the line-to-line voltage.
    Qc_var = (V_target_v**2) / Xc
    
    # Convert the result from VAR to MVAR for the final output
    Qc_mvar = Qc_var / 1e6
    
    # --- Output Results ---
    print("Reactive Power Compensation Analysis for Bus 4")
    print("-" * 50)
    print("This calculation determines the reactive power supplied by the specified shunt capacitor bank to maintain the voltage at Bus 4.")
    print("\nGiven values:")
    print(f"Target Voltage (V) = {V_target_kv:.0f} kV")
    print(f"Shunt Capacitor Reactance (Xc) = {Xc:.1f} Ω")
    
    print("\nThe formula used is: Qc = V^2 / Xc")
    print("\nCalculation with the given numbers:")
    # As requested, printing the numbers in the final equation
    print(f"Qc = ({V_target_v:.0f} V)^2 / {Xc:.1f} Ω")
    print(f"Qc = {Qc_var:,.0f} VAR")
    print(f"Qc = {Qc_mvar:,.0f} MVAR")
    print("-" * 50)
    
    print(f"\nThe reactive power needed (supplied by the capacitor bank) to maintain the bus voltage at 400 kV is {Qc_mvar:,.0f} MVAR.")

# Execute the function
calculate_reactive_power_compensation()