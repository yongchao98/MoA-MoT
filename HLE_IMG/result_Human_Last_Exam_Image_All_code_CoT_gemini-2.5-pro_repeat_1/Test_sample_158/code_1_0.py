import math

def calculate_load_voltage():
    """
    Calculates the voltage across the load resistor in the given RF rectifier circuit.
    """
    # --- Given Parameters ---
    # Input power in Watts
    P_in = 10e-3  # 10 mW
    # Load resistance in Ohms
    R_L = 2.7e3   # 2.7 kOhm
    # Characteristic impedance in Ohms
    Z_0 = 50
    # Operating frequency in Hz (determines which Q to read from the graph)
    f = 0.8e9     # 0.8 GHz = 800 MHz
    # Quality factor of the inductors (read from the graph at 800 MHz)
    Q_L = 100
    # Quality factor of the capacitors (given)
    Q_C = 150

    # --- Analysis based on the plan ---

    # Step 1: Calculate the effective input resistance of the voltage doubler rectifier.
    # R_in_rect ≈ R_L / 2
    R_in_rect = R_L / 2

    # Step 2: Calculate the loaded Q (Q_loaded) required for the impedance matching network.
    # The network matches the high impedance R_in_rect to the low impedance Z_0.
    # Q_loaded = sqrt(R_high / R_low - 1)
    Q_loaded = math.sqrt(R_in_rect / Z_0 - 1)

    # Step 3: Calculate the efficiency of the matching network (eta_match).
    # We assume the unloaded Q (Q_u) is dominated by the inductor's Q_L.
    # The efficiency is approximated as eta ≈ 1 - 2 * (Q_loaded / Q_u).
    Q_u = Q_L
    eta_match = 1 - 2 * (Q_loaded / Q_u)

    # Step 4: Calculate the DC power delivered to the load.
    # Assuming ideal diode conversion, P_DC is the input power scaled by the matching efficiency.
    P_DC = P_in * eta_match

    # Step 5: Calculate the final voltage across the load resistor R_L.
    # V_RL = sqrt(P_DC * R_L)
    V_RL = math.sqrt(P_DC * R_L)

    # --- Output the results step-by-step ---
    print("--- Calculation Steps ---")
    print(f"1. Approximated rectifier input resistance (R_in_rect = R_L / 2): {R_L:.1f} Ω / 2 = {R_in_rect:.1f} Ω")
    print(f"2. Required loaded Q (Q_loaded = sqrt(R_in_rect / Z_0 - 1)): sqrt({R_in_rect:.1f} / {Z_0:.1f} - 1) = {Q_loaded:.3f}")
    print(f"3. Matching network efficiency (eta_match = 1 - 2 * Q_loaded / Q_u): 1 - 2 * {Q_loaded:.3f} / {Q_u:.1f} = {eta_match:.4f}")
    print(f"4. DC Power to load (P_DC = P_in * eta_match): {P_in:.4f} W * {eta_match:.4f} = {P_DC:.6f} W")
    
    print("\n--- Final Voltage Calculation ---")
    print("The final voltage is calculated using the formula: V_RL = sqrt(P_DC * R_L)")
    # The final equation with all numbers plugged in
    final_equation = f"V_RL = sqrt({P_DC:.6f} W * {R_L:.1f} Ω)"
    print("Plugging in the calculated values:")
    print(final_equation)
    
    print(f"\nResult: The calculated voltage across the load resistor is {V_RL:.2f} V.")

# Execute the function
calculate_load_voltage()