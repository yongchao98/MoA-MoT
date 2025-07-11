import math

def calculate_load_voltage():
    """
    Calculates the voltage across the load resistor in the given RF harvesting circuit.

    This calculation is based on a model that estimates the total power conversion
    efficiency (PCE) by considering losses from passive components (inductors and
    capacitors) and an assumed baseline efficiency for the rectifier diodes.
    """

    # --- Step 1: Define Constants and Given Values ---
    f = 0.8e9  # Operating frequency in Hz (0.8 GHz)
    P_in = 10e-3  # Input power in Watts (10 mW)
    Z0 = 50.0  # Characteristic impedance in Ohms
    C1 = 0.4e-12  # Capacitance of C1 in Farads
    L1 = 43e-9  # Inductance of L1 in Henrys
    C2 = 0.2e-12  # Capacitance of C2 in Farads
    L2 = 39e-9  # Inductance of L2 in Henrys
    RL = 2.7e3  # Load resistance in Ohms (2.7 kOhm)

    # --- Step 2: Component Quality Factors ---
    # From the graph at 800 MHz, the Quality factor of the "Proposed" inductor is ~100.
    # Assumption: This Q-factor applies to all inductors in the circuit.
    QL = 100.0
    # The capacitor quality factor is given.
    QC = 150.0

    print("--- Assumptions ---")
    print(f"Inductor Quality Factor (QL) assumed for all inductors: {QL}")
    print(f"Capacitor Quality Factor (QC): {QC}")
    print("Assuming perfect input matching (no reflection loss).")
    print("Assuming a baseline rectifier efficiency of 65% to account for diode losses.\n")

    # --- Step 3: Calculate Unloaded Q of Sub-circuits ---
    omega = 2 * math.pi * f  # Angular frequency in rad/s

    # Series arm (L1, C1)
    X_L1 = omega * L1
    X_C1 = -1 / (omega * C1)
    R_sL1 = X_L1 / QL
    R_sC1 = abs(X_C1) / QC
    R_s_total = R_sL1 + R_sC1
    X_s_total = X_L1 + X_C1
    Q_series_arm = abs(X_s_total) / R_s_total

    # Parallel arm (L2, C2)
    X_L2 = omega * L2
    X_C2 = -1 / (omega * C2)
    # Convert series parasitic resistance to equivalent parallel resistance for loss calculation
    R_pL2 = (X_L2 / QL) * (1 + QL**2)
    R_pC2 = (abs(X_C2) / QC) * (1 + QC**2)
    R_p_loss = (R_pL2 * R_pC2) / (R_pL2 + R_pC2) # Resistors in parallel
    # Reactance of the ideal L2 || C2 tank
    X_p_tank = (X_L2 * X_C2) / (X_L2 + X_C2)
    Q_parallel_arm = R_p_loss / abs(X_p_tank)
    
    # Q of the matching inductor is simply its component Q
    Q_match_inductor = QL

    # --- Step 4: Calculate Passive Network Efficiency ---
    # Model efficiency as a product of (1 - 1/Q) for each lossy section
    eta_match_inductor = 1 - (1 / Q_match_inductor)
    eta_series_arm = 1 - (1 / Q_series_arm)
    eta_parallel_arm = 1 - (1 / Q_parallel_arm)
    
    eta_passive = eta_match_inductor * eta_series_arm * eta_parallel_arm

    # --- Step 5: Assume a Baseline Rectifier Efficiency ---
    # This accounts for diode non-idealities (e.g., forward voltage drop, junction capacitance)
    eta_base_rectifier = 0.65  # 65%

    # --- Step 6: Calculate Total Power Conversion Efficiency (PCE) ---
    total_efficiency = eta_passive * eta_base_rectifier
    
    # --- Step 7: Calculate Load Power and Voltage ---
    P_load = P_in * total_efficiency
    V_load = math.sqrt(P_load * RL)
    
    # --- Print Results Step-by-Step ---
    print("--- Calculations ---")
    print(f"Operating Frequency (f): {f/1e6:.1f} MHz")
    print(f"Angular Frequency (w): {omega/1e9:.3f} Grad/s")
    print(f"Input Power (Pin): {P_in * 1000:.1f} mW\n")
    
    print("1. Q-factor of the Series Arm (L1, C1):")
    print(f"   - Combined Resistance (Rs_tot): {R_s_total:.3f} Ohms")
    print(f"   - Combined Reactance (Xs_tot): {X_s_total:.3f} Ohms")
    print(f"   - Unloaded Q (Qs): {Q_series_arm:.3f}\n")
    
    print("2. Q-factor of the Parallel Arm (L2, C2):")
    print(f"   - Equivalent Parallel Loss Resistance (Rp_loss): {R_p_loss/1e3:.3f} kOhms")
    print(f"   - Tank Reactance (Xp_tank): {X_p_tank:.3f} Ohms")
    print(f"   - Unloaded Q (Qp): {Q_parallel_arm:.3f}\n")
    
    print("3. Efficiency Calculation:")
    print(f"   - Efficiency of Passive Network (eta_passive): {eta_passive:.4f}")
    print(f"   - Assumed Baseline Rectifier Efficiency (eta_base): {eta_base_rectifier:.4f}")
    print(f"   - Total Power Conversion Efficiency (PCE): {total_efficiency:.4f}\n")

    print("4. Final Voltage Calculation:")
    print(f"   - Power Delivered to Load (P_load): {P_load * 1000:.3f} mW")
    print(f"   - Voltage Across Load (V_load): {V_load:.3f} V\n")

    # Final equation format requested by the user
    print("--- Final Equation ---")
    print(f"V_L = sqrt(P_in * total_efficiency * R_L)")
    print(f"V_L = sqrt({P_in:.4f} W * {total_efficiency:.4f} * {RL:.1f} Ohms)")
    print(f"V_L = sqrt({P_load:.6f} * {RL:.1f})")
    print(f"V_L = {V_load:.3f} V")
    return V_load

# Execute the calculation
final_voltage = calculate_load_voltage()
# The final answer format, as requested.
# Using 2 significant figures for the final answer as is common in engineering.
print(f"<<<{final_voltage:.2f}>>>")