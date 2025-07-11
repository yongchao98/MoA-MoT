import math

def calculate_reactive_power():
    """
    Calculates the reactive power supplied by a shunt capacitor bank.

    The problem asks to determine the reactive power needed to maintain the bus voltage
    at 400 kV, specifying that the compensation is provided by a shunt capacitor
    bank with a reactance of Xc = 5 Ω.

    The reactive power (Qc) supplied by a three-phase shunt capacitor bank is
    calculated using the formula:
    Qc = V_LL^2 / Xc
    where V_LL is the line-to-line voltage and Xc is the capacitor's reactance.

    Other information in the problem, such as the 5% voltage drop, is not used as
    there is insufficient data to perform a full power flow analysis. The calculation
    is based on the explicit values provided for the capacitor and the target voltage.
    """

    # Given parameters
    V_ll_kv = 400  # Line-to-line voltage in kV
    Xc = 5         # Capacitor reactance in Ohms

    # Convert voltage to Volts for calculation
    V_ll_volts = V_ll_kv * 1000

    # Calculate reactive power in VAR
    Qc_var = V_ll_volts**2 / Xc

    # Convert result to MVAR
    Qc_mvar = Qc_var / 1_000_000

    print("Step 1: Identify the formula and given values.")
    print("Formula for reactive power from a shunt capacitor: Qc = V_LL^2 / Xc")
    print(f"Given Line-to-Line Voltage (V_LL): {V_ll_kv} kV = {int(V_ll_volts)} V")
    print(f"Given Capacitor Reactance (Xc): {Xc} Ω\n")

    print("Step 2: Perform the calculation.")
    # Printing the equation with all the numbers as requested
    print(f"Qc = ({int(V_ll_volts)})^2 / {Xc}")
    print(f"Qc = {int(V_ll_volts**2)} / {Xc}")
    print(f"Qc = {int(Qc_var)} VAR\n")


    print("Step 3: Convert the result to MegaVAR (MVAR).")
    print(f"Qc = {int(Qc_var)} / 1,000,000")
    print(f"Qc = {int(Qc_mvar)} MVAR\n")

    print("The reactive power needed to be supplied by the capacitor bank is {:.0f} MVAR.".format(Qc_mvar))

if __name__ == "__main__":
    calculate_reactive_power()