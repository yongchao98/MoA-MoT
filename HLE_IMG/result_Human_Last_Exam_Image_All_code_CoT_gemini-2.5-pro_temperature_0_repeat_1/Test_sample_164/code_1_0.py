import math

def calculate_reactive_power_compensation():
    """
    Performs a reactive power compensation analysis to maintain voltage stability.

    This function calculates the reactive power needed at Bus 4 to compensate for a
    voltage drop on the line from Bus 5, maintaining the Bus 4 voltage at 400 kV.
    """
    # Given parameters
    V_target_kV = 400.0  # Target voltage at Bus 4 in kV
    voltage_drop_percent = 5.0  # Voltage drop on the line from Bus 5 to Bus 4
    X_c = 5.0  # Reactance of the shunt capacitor bank in Ohms

    # Convert target voltage to Volts
    V_target_V = V_target_kV * 1000

    # Calculate the voltage drop magnitude in Volts
    delta_V_V = (voltage_drop_percent / 100.0) * V_target_V

    # The approximate voltage drop formula for a reactive line is: Delta_V ≈ (Q * X) / V
    # The reactive power needed for compensation (Qc) must offset the reactive power
    # consumed by the line. We can find this required Qc by rearranging the formula:
    # Qc = (Delta_V * V) / X
    # We assume the effective reactance X is the given capacitor reactance X_c.
    Qc_var = (delta_V_V * V_target_V) / X_c

    # Convert reactive power from VAR to MVAR
    Qc_mvar = Qc_var / 1e6

    # Print the detailed analysis and the final equation with values
    print("--- Reactive Power Compensation Analysis ---")
    print(f"Target voltage at Bus 4 (V): {V_target_kV} kV")
    print(f"Voltage drop to compensate (ΔV): {voltage_drop_percent}% of {V_target_kV} kV = {delta_V_V / 1000} kV")
    print(f"Effective system reactance (X): {X_c} Ω")
    print("\nFormula to calculate required reactive power (Qc):")
    print("Qc = (ΔV * V) / X")
    print("\nFinal Calculation:")
    # Using integers for cleaner printing of the equation
    print(f"Qc = ({int(delta_V_V)} V * {int(V_target_V)} V) / {X_c} Ω")
    print(f"Qc = {int(delta_V_V * V_target_V)} / {X_c} VAR")
    print(f"Qc = {int(Qc_var)} VAR")
    print(f"Qc = {int(Qc_mvar)} MVAR")

# Execute the function
calculate_reactive_power_compensation()

# The final numerical answer
final_answer = ( (voltage_drop_percent / 100.0) * (V_target_kV * 1000) * (V_target_kV * 1000) ) / (X_c * 1e6)
print(f"\n<<<1600>>>")