import math

def calculate_load_voltage():
    """
    Calculates the voltage across the load resistor in the given RF circuit.
    """
    # Step 1: Extract Parameters
    P_in = 10e-3  # Input Power in Watts (10 mW)
    f = 0.8e9     # Operating Frequency in Hz (0.8 GHz)
    R_L = 2.7e3   # Load Resistance in Ohms (2.7 kOhms)
    L1 = 43e-9    # Inductance L1 in Henries (43 nH)
    C1 = 0.4e-12  # Capacitance C1 in Farads (0.4 pF)
    L2 = 39e-9    # Inductance L2 in Henries (39 nH)
    C2 = 0.2e-12  # Capacitance C2 in Farads (0.2 pF)
    Q_C = 150     # Quality factor of capacitors
    # From the graph at 800 MHz, the Quality Factor for the inductor is at its peak.
    Q_L = 100     # Quality factor of inductors

    print("--- Step 1: Initial Parameters ---")
    print(f"Input Power (Pin): {P_in * 1000} mW")
    print(f"Frequency (f): {f / 1e6} MHz")
    print(f"Load Resistance (RL): {R_L / 1000} kOhms")
    print(f"Inductor Quality Factor (QL): {Q_L}")
    print(f"Capacitor Quality Factor (QC): {Q_C}\n")

    # Step 2: Calculate Angular Frequency and Effective Load Resistance
    omega = 2 * math.pi * f
    # For an ideal voltage doubler, the effective RF input resistance is RL/2
    R_in_load = R_L / 2

    print("--- Step 2: Intermediate Calculations ---")
    print(f"Angular Frequency (w): {omega:.4g} rad/s")
    print(f"Effective RF Load Resistance (R_in_load): {R_in_load} Ohms\n")

    # Step 3: Calculate Equivalent Series Resistance (ESR) for each component
    Rs_L1 = (omega * L1) / Q_L
    Rs_C1 = 1 / (omega * C1 * Q_C)
    Rs_L2 = (omega * L2) / Q_L
    Rs_C2 = 1 / (omega * C2 * Q_C)

    print("--- Step 3: Component Loss Calculation (ESR) ---")
    print(f"ESR of L1 (Rs_L1): {Rs_L1:.4f} Ohms")
    print(f"ESR of C1 (Rs_C1): {Rs_C1:.4f} Ohms")
    print(f"ESR of L2 (Rs_L2): {Rs_L2:.4f} Ohms")
    print(f"ESR of C2 (Rs_C2): {Rs_C2:.4f} Ohms\n")

    # Step 4: Calculate Total Loss Resistance
    R_loss = Rs_L1 + Rs_C1 + Rs_L2 + Rs_C2
    print("--- Step 4: Total Loss Resistance ---")
    print(f"Total Loss Resistance (R_loss): {R_loss:.4f} Ohms\n")

    # Step 5: Calculate Conduction Efficiency
    # Assuming perfect matching (eta_M=1), so total efficiency = conduction efficiency
    eta_total = R_in_load / (R_in_load + R_loss)
    print("--- Step 5: Efficiency Calculation ---")
    print(f"Conduction Efficiency (eta_c): {eta_total:.4f}\n")
    
    # Step 6: Calculate Output Power
    P_out = P_in * eta_total
    print("--- Step 6: Output Power Calculation ---")
    print(f"Output Power (Pout): {P_out * 1000:.4f} mW\n")

    # Step 7: Calculate Load Voltage
    V_DC = math.sqrt(P_out * R_L)
    print("--- Step 7: Final Voltage Calculation ---")
    print("The final equation for the load voltage V_DC is:")
    print(f"V_DC = sqrt(P_in * eta_total * R_L)")
    print(f"V_DC = sqrt({P_in:.3f} W * {eta_total:.4f} * {R_L:.0f} Ohms)")
    print(f"V_DC = sqrt({P_out:.6f} W * {R_L:.0f} Ohms)")
    print(f"V_DC = sqrt({P_out * R_L:.4f} V^2)")
    print(f"Calculated V_DC = {V_DC:.4f} V")
    
    return V_DC

# Run the calculation and print the final result
final_voltage = calculate_load_voltage()
print(f"\nThe final calculated voltage across the load resistor is {final_voltage:.3f} V.")
print(f"<<<{final_voltage:.3f}>>>")