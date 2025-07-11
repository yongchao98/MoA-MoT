import math
import cmath

def solve_circuit_voltage():
    """
    Calculates the total voltage across the load capacitor in the given circuit,
    considering harmonics and parasitic effects.
    """

    # 1. Define given and assumed parameters
    V_RF = 1.0  # V (Fundamental amplitude)
    f0 = 915e6  # Hz (Fundamental frequency)
    R0 = 50.0   # Ohm (Base parasitic resistance)
    C_parasitic = 2e-15 # F (Parasitic capacitance)
    R_L = 8000.0 # Ohm (Load resistance)
    
    # Assumptions for the rectifier model
    # Effective input resistance of the rectifier
    R_in_rect = R_L / 2.0  
    # Power conversion efficiency (PCE)
    eta = 1.0 

    harmonics = [1, 3, 5, 7]
    total_power_in = 0.0

    print("### System Parameters and Assumptions ###")
    print(f"Fundamental Voltage Amplitude V_RF = {V_RF} V")
    print(f"Fundamental Frequency f0 = {f0/1e6} MHz")
    print(f"Load Resistance R_L = {R_L} Ohm")
    print(f"Parasitic Resistance R0 = {R0} Ohm at f0")
    print(f"Parasitic Capacitance C_parasitic = {C_parasitic*1e15} fF")
    print(f"Assumed Rectifier Input Resistance R_in_rect = {R_in_rect} Ohm")
    print(f"Assumed Power Conversion Efficiency eta = {eta*100}%\n")

    # 2. Calculate power contribution from each harmonic
    print("### Calculating Power Delivered by Each Harmonic ###")
    for n in harmonics:
        # Calculate harmonic-specific values
        V_n = V_RF * (0.9**((n - 1) // 2))
        f_n = n * f0
        w_n = 2 * math.pi * f_n
        R_p_n = R0 * (n**2)

        # Calculate impedance of the effective load (R_in_rect || C_parasitic) using complex numbers
        # Z_load_eff = 1 / (1/R_in_rect + j*w_n*C_parasitic)
        Z_load_eff = 1 / (1/R_in_rect + 1j * w_n * C_parasitic)
        
        # Calculate total impedance of the signal path
        Z_total = R_p_n + Z_load_eff
        
        # Calculate voltage across the effective load using the voltage divider rule
        V_load_eff = V_n * (Z_load_eff / Z_total)
        
        # Power delivered is based on the voltage across the resistive part (R_in_rect)
        # P = V_peak^2 / (2*R). The voltage across R_in_rect is V_load_eff.
        power_n = (abs(V_load_eff)**2) / (2 * R_in_rect)
        total_power_in += power_n

        print(f"--- Harmonic n={n} ---")
        print(f"  Voltage Amplitude V_{n} = {V_n:.3f} V")
        print(f"  Parasitic Resistance R_p({n}) = {R_p_n:.1f} Ohm")
        print(f"  Voltage at Rectifier Input |V_load_eff| = {abs(V_load_eff):.4f} V")
        print(f"  Power Delivered P_in({n}) = {power_n * 1e6:.2f} uW")

    # 3. Calculate final DC output voltage
    V_DC_squared = eta * total_power_in * R_L
    V_DC = math.sqrt(V_DC_squared)

    print("\n### Final Voltage Calculation ###")
    print(f"Total AC power delivered to rectifier: P_in_total = {total_power_in:.6f} W")
    print("The final DC voltage is calculated using: V_DC = sqrt(eta * P_in_total * R_L)")
    print(f"V_DC = sqrt({eta} * {total_power_in:.6f} * {R_L})")
    print(f"V_DC = sqrt({V_DC_squared:.4f})")
    print(f"Final calculated voltage across C_L, V_DC = {V_DC:.3f} V")

solve_circuit_voltage()
<<<1.486>>>