import math

def calculate_load_voltage():
    """
    Calculates the voltage across the load resistor in the given RF circuit.
    """
    # --- Step 1: Define Given Parameters ---
    # Input power in Watts (10 mW)
    P_in = 10e-3
    # Load resistance in Ohms (2.7 kOhm)
    R_L = 2.7e3
    # Operating frequency in Hz (0.8 GHz = 800 MHz)
    f = 0.8e9
    # From Figure (b), the quality factor of the inductor at 800 MHz is read from the graph.
    # The red dashed line (Quality factor) at 800 MHz on the x-axis corresponds to 90 on the right y-axis.
    Q_L = 90
    # The quality factor of the capacitor is given in the text.
    Q_C = 150

    print("--- Analysis of Given Information ---")
    print(f"Input Power (P_in): {P_in * 1000} mW")
    print(f"Load Resistance (R_L): {R_L / 1000} kOhm")
    print(f"Operating Frequency (f): {f / 1e6} MHz")
    print(f"Inductor Quality Factor (Q_L) from graph: {Q_L}")
    print(f"Capacitor Quality Factor (Q_C): {Q_C}")
    print("-" * 35 + "\n")

    # --- Step 2: Determine the Power Transfer Efficiency (eta) ---
    # The statement "tuned for optimal power transfer" suggests the circuit is critically coupled.
    # Under critical coupling, the power transfer efficiency through a passive resonant network is 50%,
    # as half the power is dissipated in the components and half is delivered to the load.
    # We assume reflection losses are minimized to zero by the matching circuit, and that the diode is ideal,
    # as the problem focuses on losses due to component Q-factors.
    efficiency = 0.5
    
    print("--- Efficiency Calculation ---")
    print("Assumption: 'Optimal power transfer' implies the system is critically coupled.")
    print("The efficiency (eta) of a critically coupled passive network is 50%.")
    print(f"Overall Efficiency (eta) = {efficiency:.2f}\n")

    # --- Step 3: Calculate the Output DC Power (P_out) ---
    P_out = P_in * efficiency

    print("--- Output Power Calculation ---")
    print(f"The output DC power (P_out) is calculated as: P_out = P_in * eta")
    print(f"P_out = {P_in:.3f} W * {efficiency:.2f}")
    print(f"P_out = {P_out:.4f} W\n")

    # --- Step 4: Calculate the Voltage across the Load (V_L) ---
    # The DC power dissipated by the load is P_out = V_L^2 / R_L
    # Rearranging for V_L gives V_L = sqrt(P_out * R_L)
    V_L = math.sqrt(P_out * R_L)
    intermediate_value = P_out * R_L

    print("--- Output Voltage Calculation ---")
    print(f"The voltage across the load (V_L) is calculated as: V_L = sqrt(P_out * R_L)")
    print(f"V_L = sqrt({P_out:.4f} W * {R_L:.0f} Ohms)")
    print(f"V_L = sqrt({intermediate_value:.2f}) V")
    print(f"V_L = {V_L:.3f} V")

    # --- Final Answer ---
    print("\n<<<3.674>>>")

calculate_load_voltage()