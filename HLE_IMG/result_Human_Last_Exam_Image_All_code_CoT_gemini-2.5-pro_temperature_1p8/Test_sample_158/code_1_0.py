import math

def calculate_load_voltage():
    """
    Calculates the DC voltage across the load resistor in an RF energy harvesting circuit.
    The calculation follows a detailed model including matching network efficiency,
    resonant tank losses, and diode forward voltage drop.
    """
    # Step 1: Define constants and parameters
    f = 0.8e9  # Operating frequency in Hz (0.8 GHz)
    P_in = 0.01  # Input power in Watts (10 mW)
    RL = 2700.0  # Load resistance in Ohms (2.7 kOhm)
    Z0 = 50.0  # Characteristic impedance in Ohms
    Q_L = 98.0  # Inductor quality factor from graph at 800 MHz
    Q_C = 150.0  # Capacitor quality factor
    L1 = 43.0e-9  # Inductance L1 in Henrys
    C1 = 0.4e-12  # Capacitance C1 in Farads
    L2 = 39.0e-9  # Inductance L2 in Henrys
    C2 = 0.2e-12  # Capacitance C2 in Farads
    Vf = 0.3  # Estimated diode forward voltage in Volts

    omega = 2 * math.pi * f  # Angular frequency in rad/s

    print(f"--- Step 1: Initial Parameters ---")
    print(f"Operating Frequency (f): {f/1e9:.2f} GHz")
    print(f"Input Power (Pin): {P_in*1000:.1f} mW")
    print(f"Inductor Q-factor (QL): {Q_L}")
    print(f"Capacitor Q-factor (QC): {Q_C}\n")

    # Step 2: Model rectifier input impedance
    n_multiplier = 2  # Voltage doubler
    R_rect = RL / (n_multiplier**2)
    print(f"--- Step 2: Rectifier Input Impedance ---")
    print(f"Estimated rectifier input resistance (R_rect): {R_rect:.2f} Ohms\n")

    # Step 3: Calculate matching network efficiency
    # Loaded Q for the matching network transforming R_rect to Z0
    Q_load = math.sqrt(R_rect / Z0 - 1)
    # Efficiency of the matching network dominated by inductor Q
    eta_match = 1 - Q_load / Q_L
    P_rect_in = P_in * eta_match
    print(f"--- Step 3: Matching Network ---")
    print(f"Matching network loaded Q (Q_load): {Q_load:.4f}")
    print(f"Matching network efficiency (eta_match): {eta_match:.4f}")
    print(f"Power at rectifier input (P_rect_in): {P_rect_in*1000:.4f} mW\n")

    # Step 4: Calculate power loss in resonant tanks
    V_rms_rect = math.sqrt(P_rect_in * R_rect)
    I_rms_rect = V_rms_rect / R_rect
    
    # Loss in series tank L1-C1
    R_sL1 = (omega * L1) / Q_L
    X_C1 = 1 / (omega * C1)
    R_sC1 = X_C1 / Q_C
    R_tank1_series = R_sL1 + R_sC1
    P_loss_T1 = (I_rms_rect**2) * R_tank1_series
    
    # Loss in parallel tank L2-C2
    R_pL2 = (omega * L2) * Q_L
    X_C2 = 1 / (omega * C2)
    R_pC2 = X_C2 * Q_C
    R_tank2_parallel = (R_pL2 * R_pC2) / (R_pL2 + R_pC2)
    P_loss_T2 = (V_rms_rect**2) / R_tank2_parallel

    print(f"--- Step 4: Resonant Tank Losses ---")
    print(f"RMS Voltage at rectifier (V_rms_rect): {V_rms_rect:.4f} V")
    print(f"RMS Current at rectifier (I_rms_rect): {I_rms_rect*1000:.4f} mA")
    print(f"Series Tank 1 parasitic resistance: {R_tank1_series:.4f} Ohms")
    print(f"Power loss in Tank 1 (P_loss_T1): {P_loss_T1*1000:.4f} mW")
    print(f"Parallel Tank 2 parasitic resistance: {R_tank2_parallel/1000:.4f} kOhms")
    print(f"Power loss in Tank 2 (P_loss_T2): {P_loss_T2*1000:.4f} mW\n")

    # Step 5: Power available for diode conversion
    P_prime_rect = P_rect_in - P_loss_T1 - P_loss_T2
    print(f"--- Step 5: Power for Rectification ---")
    print(f"Power delivered to ideal diodes + load (P'_rect): {P_prime_rect*1000:.4f} mW\n")

    # Step 6: Solve for load voltage VL
    # Equation: a*VL^2 + b*VL + c = 0
    a = 1.0
    b = 2 * Vf
    c = -P_prime_rect * RL
    
    # Using quadratic formula: VL = (-b + sqrt(b^2 - 4ac)) / 2a
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        print("Error: No real solution for VL.")
        return
        
    VL = (-b + math.sqrt(discriminant)) / (2 * a)

    print(f"--- Step 6: Final Voltage Calculation ---")
    print(f"Solving the quadratic equation for VL: a*VL^2 + b*VL + c = 0")
    print(f"Where:")
    print(f"a = {a}")
    print(f"b = 2 * Vf = {b:.4f}")
    print(f"c = -P'_rect * RL = {c:.4f}")
    print(f"\nThe calculated voltage across the load RL is: {VL:.4f} V")

    return VL

# Run the calculation and store the final answer
final_voltage = calculate_load_voltage()
# We are asked to provide the final numerical answer in the specified format
# The value is ~4.69
print(f"\nFinal Answer: {final_voltage:.2f} V")
print("<<<4.69>>>")
