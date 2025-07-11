import cmath
import math

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # --- 1. Define Constants and Assumptions ---
    # From problem description and image table
    V_RF = 1.0          # Peak voltage of the fundamental frequency [V]
    f0 = 915e6          # Fundamental frequency [Hz]
    R0 = 50.0           # Base parasitic resistance at f0 [Ohm]
    C_parasitic = 2e-15 # Parasitic capacitance [F]
    
    # Assumption: The rectifier input impedance is matched to the system.
    # The R0 = 50 Ohm parameter strongly suggests a 50 Ohm system.
    Z_rectifier = 50.0 + 0j # Rectifier input impedance [Ohm]

    print("--- System Model and Assumptions ---")
    print(f"Fundamental Frequency (f0): {f0 / 1e6} MHz")
    print(f"Parasitic Resistance R_parasitic(f) = {R0} * (f/f0)^2 Ohms")
    print(f"Parasitic Capacitance (C_parasitic): {C_parasitic * 1e15} fF")
    print(f"Assumed Rectifier Input Impedance (Z_rectifier): {Z_rectifier.real} Ohms")
    print("Efficiency Definition: P_useful / P_in_total")
    print("  - P_useful: Power at f0 delivered to Z_rectifier.")
    print("  - P_in_total: Total power from the source across all harmonics.\n")

    # --- 2. Calculate Harmonic Voltages ---
    harmonics = [1, 3, 5, 7]
    V = {}
    V[1] = V_RF
    # Voltage drops by 10% for each higher harmonic
    for i in range(len(harmonics) - 1):
        V[harmonics[i+1]] = V[harmonics[i]] * 0.9

    # --- 3. Loop Through Harmonics to Calculate Power ---
    P_in_total = 0.0
    P_rect_1 = 0.0
    
    # To store power values for final equation output
    power_sources = {}

    print("--- Power Calculation per Harmonic ---")
    for n in harmonics:
        V_n = V[n]
        f_n = n * f0
        omega_n = 2 * math.pi * f_n
        
        # Parasitic resistance for the nth harmonic
        R_p_n = R0 * (n**2)
        
        # Impedance of the shunt parasitic capacitor
        Z_Cp_n = 1 / (1j * omega_n * C_parasitic)
        
        # Impedance of the parallel combination (rectifier || C_parasitic)
        Z_L_prime_n = 1 / (1/Z_rectifier + 1/Z_Cp_n)
        
        # Total impedance seen by the source
        Z_total_n = R_p_n + Z_L_prime_n
        
        # Average power drawn from the source for this harmonic
        # P = |V_rms|^2 * Re(Z) / |Z|^2, where V_rms = V_peak / sqrt(2)
        P_source_n = (V_n**2 / 2) * (Z_total_n.real / (abs(Z_total_n)**2))
        power_sources[n] = P_source_n
        P_in_total += P_source_n
        
        # If this is the fundamental, calculate the "useful" power delivered to the rectifier
        if n == 1:
            # Power delivered to Z_L_prime_n is dissipated in its real part
            # Current I_n_rms = (V_n/sqrt(2)) / Z_total_n
            # P = |I_n_rms|^2 * Re(Z_L_prime_n)
            P_rect_1 = (V_n**2 / 2) * (Z_L_prime_n.real / (abs(Z_total_n)**2))

        print(f"Harmonic n={n}:")
        print(f"  V_{n} = {V_n:.3f} V")
        print(f"  P_source_{n} = {P_source_n * 1000:.4f} mW")

    # --- 4. Calculate Final Efficiency ---
    efficiency = P_rect_1 / P_in_total

    print("\n--- Overall System Efficiency Calculation ---")
    print(f"Useful Power (P_rect_1): {P_rect_1 * 1000:.4f} mW")
    print(f"Total Input Power (P_in_total): {P_in_total * 1000:.4f} mW")
    
    print("\nFinal Equation:")
    p_s1 = power_sources[1]
    p_s3 = power_sources[3]
    p_s5 = power_sources[5]
    p_s7 = power_sources[7]
    
    print(f"Efficiency = P_rect_1 / (P_source_1 + P_source_3 + P_source_5 + P_source_7)")
    print(f"Efficiency = {P_rect_1:.6f} / ({p_s1:.6f} + {p_s3:.6f} + {p_s5:.6f} + {p_s7:.6f})")
    print(f"Efficiency = {P_rect_1:.6f} / {P_in_total:.6f}")
    print(f"Efficiency = {efficiency:.4f}")
    print(f"Efficiency = {efficiency * 100:.2f} %")
    
    return efficiency

# Run the calculation and print the result
final_efficiency = calculate_efficiency()
print(f"\n<<<>>>\n") # Suppress the final raw number return and add the answer tag as requested
# The logic above already printed the answer.
# To satisfy the format, let's print just the number.
print(f"<<<{final_efficiency:.4f}>>>")

calculate_efficiency()