import math

def calculate_load_voltage():
    """
    Calculates the DC voltage across the load resistor in the given RF rectifier circuit.
    
    The calculation follows a model where the circuit's overall efficiency is determined by
    the quality factors of its components and the required impedance transformation.
    """
    # Step 1: Extract Given Parameters
    Pin = 10e-3  # Input Power in Watts (10 mW)
    RL = 2.7e3   # Load Resistance in Ohms (2.7 KΩ)
    Z0 = 50      # Characteristic Impedance in Ohms
    f = 0.8e9    # Operating Frequency in Hz (0.8 GHz)
    Q_C = 150    # Quality factor of capacitors

    # Step 2: Determine Inductor Quality Factor from the graph
    # From Figure (b), at 800 MHz, the Quality Factor (dashed line) is approximately 90.
    Q_L = 90
    
    print("--- Input Parameters ---")
    print(f"Input Power (Pin): {Pin} W")
    print(f"Load Resistance (RL): {RL} Ohms")
    print(f"Source Impedance (Z0): {Z0} Ohms")
    print(f"Inductor Quality Factor (Q_L): {Q_L}")
    print(f"Capacitor Quality Factor (Q_C): {Q_C}\n")

    # Step 3: Model the circuit for efficiency calculation
    # Estimate the rectifier's effective AC input resistance
    R_rect = RL / 2
    print("--- Efficiency Calculation ---")
    print(f"Step 3a: Estimated Rectifier AC Input Resistance (R_rect = RL / 2): {R_rect:.2f} Ohms")

    # Step 4: Calculate Loaded and Unloaded Q-Factors
    # Calculate the required loaded Q (Q_L_loaded) for impedance matching from Z0 to R_rect
    Q_L_loaded = math.sqrt(R_rect / Z0 - 1)
    print(f"Step 4a: Calculated Loaded Q (Q_L_loaded): {Q_L_loaded:.4f}")

    # Calculate the combined unloaded Q (Q_U) of the passive network.
    # The network includes 3 inductors (1 for matching, 2 for rectifier circuit)
    # and 2 capacitors (for rectifier circuit).
    # 1/Q_U = (1/Q_L + 1/Q_L + 1/Q_L) + (1/Q_C + 1/Q_C)
    inv_Q_U = 3 / Q_L + 2 / Q_C
    Q_U = 1 / inv_Q_U
    print(f"Step 4b: Calculated Unloaded Q (Q_U): {Q_U:.4f}")

    # Step 5: Calculate Overall Power Conversion Efficiency (eta)
    eta = 1 - (Q_L_loaded / Q_U)
    print(f"Step 5: Overall Power Conversion Efficiency (eta = 1 - Q_L_loaded / Q_U): {eta:.4f}\n")

    # Step 6: Calculate DC Output Power (Pdc)
    P_dc = Pin * eta

    # Step 7: Calculate DC Voltage across the Load (Vdc)
    V_dc = math.sqrt(P_dc * RL)
    
    print("--- Final Calculation ---")
    print(f"The final voltage V_dc is calculated using the formula: V_dc = sqrt(Pin * eta * RL)")
    print(f"V_dc = sqrt({Pin:.3f} W * {eta:.4f} * {RL:.1f} Ohms)")
    print(f"V_dc = sqrt({P_dc:.6f} W * {RL:.1f} Ohms)")
    print(f"V_dc = {V_dc:.4f} V")
    
    # Returning the final answer rounded to two decimal places
    return round(V_dc, 2)

final_voltage = calculate_load_voltage()
print(f"\nThe calculated voltage across the load RL is approximately {final_voltage} V.")
print(f"<<<{final_voltage}>>>")
