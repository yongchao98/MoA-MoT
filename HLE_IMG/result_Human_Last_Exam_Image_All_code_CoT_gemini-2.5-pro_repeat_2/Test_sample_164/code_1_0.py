import numpy as np

def calculate_reactive_power_compensation():
    """
    Calculates the reactive power supplied by a shunt capacitor bank.

    The problem asks for the reactive power needed to maintain the voltage at Bus 4
    at 400 kV, using a shunt capacitor with a given reactance. While the problem
    provides extensive system details, key data for a full power flow analysis
    (like load values and complete line impedances) are missing.

    Therefore, the most direct solution is to calculate the reactive power
    supplied by the specified capacitor bank at the target voltage.

    The formula for reactive power (Qc) from a shunt capacitor is:
    Qc = V^2 / Xc
    """

    # Given values
    V_ll_kV = 400  # Line-to-line voltage at Bus 4 in kV
    Xc = 5         # Reactance of the shunt capacitor bank in Ohms

    # Convert voltage from kV to V
    V_ll_V = V_ll_kV * 1000

    # Calculate reactive power in VAR
    # This formula gives the total three-phase reactive power
    Qc_VAR = V_ll_V**2 / Xc

    # Convert reactive power from VAR to MVAR for standard representation
    Qc_MVAR = Qc_VAR / 1e6

    # Print the explanation and the final equation
    print("Reactive Power Compensation Analysis for Bus 4:")
    print("-" * 50)
    print(f"Target Voltage at Bus 4 (V): {V_ll_kV} kV")
    print(f"Shunt Capacitor Reactance (Xc): {Xc} Ω")
    print("\nCalculation using the formula Qc = V^2 / Xc:")
    # Using numpy.format_float_scientific to avoid scientific notation in the large number
    V_squared_str = np.format_float_positional(V_ll_V**2, trim='-')
    print(f"Qc = ({V_ll_kV*1000:,} V)^2 / {Xc} Ω")
    print(f"Qc = {V_squared_str} V^2 / {Xc} Ω")
    Qc_VAR_str = np.format_float_positional(Qc_VAR, trim='-')
    print(f"Qc = {Qc_VAR_str:,} VAR")
    print("-" * 50)
    print(f"The reactive power supplied by the capacitor bank is {int(Qc_MVAR):,} MVAR.")
    print("-" * 50)


# Run the calculation
calculate_reactive_power_compensation()