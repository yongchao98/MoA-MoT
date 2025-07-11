import math

def calculate_reactive_power():
    """
    Calculates the reactive power supplied by a shunt capacitor bank.

    The problem asks for the reactive power needed to maintain the voltage at Bus 4
    at 400 kV, using a shunt capacitor bank with a given reactance. Despite the
    complex diagram, the problem simplifies to calculating the output of this
    specific capacitor at the target voltage.

    The formula for three-phase reactive power from a shunt capacitor is:
    Q = V_LL^2 / X_C
    where V_LL is the line-to-line voltage and X_C is the per-phase reactance.

    A convenient unit identity is used: Q (MVAR) = (V (kV))^2 / X (Ω).
    """

    # --- Given Parameters ---
    V_bus_kV = 400  # Target line-to-line voltage at Bus 4 in kV
    Xc_ohm = 5      # Reactance of the shunt capacitor bank in Ω

    # --- Calculation ---
    # Using the identity Q[MVAR] = (V[kV])^2 / X[Ω]
    Q_MVAR = (V_bus_kV**2) / Xc_ohm

    # --- Output Results ---
    print("To determine the required reactive power, we calculate the output of the specified capacitor bank at the target voltage.")
    print("The formula is Q_C = V_LL^2 / X_C.")
    print("\nGiven values:")
    print(f"  - Target Bus Voltage (V_LL): {V_bus_kV} kV")
    print(f"  - Capacitor Reactance (X_C): {Xc_ohm} Ω")
    
    print("\nFinal Equation:")
    # The prompt requires showing each number in the final equation.
    print(f"Reactive Power (MVAR) = ({V_bus_kV})^2 / {Xc_ohm}")
    print(f"Reactive Power (MVAR) = {V_bus_kV**2} / {Xc_ohm}")
    print(f"Reactive Power = {Q_MVAR} MVAR")

if __name__ == "__main__":
    calculate_reactive_power()
    # The final numerical answer extracted from the calculation.
    # Q_MVAR = (400**2) / 5 = 160000 / 5 = 32000.0
    print("\n<<<32000.0>>>")
