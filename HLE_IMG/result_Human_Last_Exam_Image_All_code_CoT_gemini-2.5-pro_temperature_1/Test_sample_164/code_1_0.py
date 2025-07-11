import math

def calculate_reactive_power():
    """
    Calculates the reactive power supplied by a shunt capacitor bank.

    The problem asks for the reactive power needed to maintain the voltage at Bus 4
    at 400 kV, using a shunt capacitor bank with a given reactance. This is calculated
    using the formula Qc = V^2 / Xc.
    """
    # --- Given Parameters ---
    # Line-to-line voltage at Bus 4 in kV
    V_LL_kv = 400.0
    # Reactance of the shunt capacitor bank in Ohms
    Xc_ohms = 5.0

    # --- Calculation ---
    # Convert voltage from kV to Volts for SI unit consistency
    V_LL_volts = V_LL_kv * 1000

    # Calculate the three-phase reactive power in VAR using the formula: Qc = (V_LL^2) / Xc
    # This formula is standard for a wye-connected bank with per-phase reactance Xc.
    Qc_var = V_LL_volts**2 / Xc_ohms

    # Convert the result from VAR to MVAR (Mega-VAR)
    Qc_mvar = Qc_var / 1e6

    # --- Output Results ---
    print("Performing reactive power compensation analysis for Bus 4.")
    print("The amount of reactive power supplied by the shunt capacitor is calculated to maintain the bus voltage at the desired level.\n")
    print("Formula: Qc = (V_LL^2) / Xc\n")
    print("Where:")
    print(f"  V_LL (Line-to-Line Voltage) = {V_LL_kv} kV")
    print(f"  Xc (Capacitor Reactance) = {Xc_ohms} Ω\n")

    print("Step-by-step calculation:")
    # Print the equation with the numbers plugged in
    print(f"1. Substitute the values into the formula:")
    print(f"   Qc = ({int(V_LL_volts)} V)^2 / {Xc_ohms} Ω")

    # Print the intermediate step
    print(f"2. Calculate the square of the voltage:")
    print(f"   Qc = {V_LL_volts**2:.2e} V^2 / {Xc_ohms} Ω")

    # Print the final result in VAR and MVAR
    print(f"3. Compute the final reactive power:")
    print(f"   Qc = {Qc_var:.2e} VAR")
    print(f"   Qc = {int(Qc_mvar)} MVAR\n")

    print(f"Therefore, the reactive power needed from the capacitor bank to maintain the voltage at {int(V_LL_kv)} kV is {int(Qc_mvar)} MVAR.")

# Execute the function to perform the analysis and print the results.
calculate_reactive_power()

# The final answer is the value calculated for the reactive power in MVAR.
final_answer = 32000
#<<<32000>>>