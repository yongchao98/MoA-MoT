import math

def calculate_reactive_power_compensation():
    """
    Calculates the reactive power needed at Bus 4 to maintain voltage stability.

    The problem asks for the reactive power needed to maintain the voltage at Bus 4
    at 400 kV, using a shunt capacitor bank with a given reactance.

    The reactive power (Q_C) supplied by a shunt capacitor is calculated using the formula:
    Q_C = V^2 / X_C
    where:
    V = Voltage at the bus
    X_C = Reactance of the capacitor bank
    """

    # Given values
    V_bus4_kV = 400  # Target voltage at Bus 4 in kV
    X_C = 5          # Reactance of the shunt capacitor in Ohms

    # Convert voltage from kV to V
    V_bus4_V = V_bus4_kV * 1000

    # Calculate the reactive power in VAR
    Q_C_VAR = V_bus4_V**2 / X_C

    # Convert reactive power from VAR to MVAR for standard representation
    Q_C_MVAR = Q_C_VAR / 1_000_000

    # Print the explanation and the calculation steps
    print("To determine the reactive power needed, we use the formula for a shunt capacitor: Q_C = V^2 / X_C")
    print("\nGiven values:")
    print(f"Target Voltage (V) at Bus 4 = {V_bus4_kV} kV")
    print(f"Capacitor Reactance (X_C) = {X_C} Ω")
    
    print("\nCalculation:")
    # Using f-string to format the equation with numbers
    print(f"Q_C = ({int(V_bus4_V)} V)^2 / {X_C} Ω")
    
    # Show intermediate step
    intermediate_result = V_bus4_V**2
    print(f"Q_C = {intermediate_result:.0f} V^2 / {X_C} Ω")

    # Show final result
    print(f"\nThe required reactive power is {int(Q_C_MVAR)} MVAR.")

calculate_reactive_power_compensation()
<<<32000>>>