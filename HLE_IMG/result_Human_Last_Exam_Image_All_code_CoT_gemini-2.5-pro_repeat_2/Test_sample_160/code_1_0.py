import cmath
import math

def calculate_output_voltage():
    """
    Calculates the total voltage across the load capacitor CL based on the provided circuit parameters and effects.
    """
    # --- 1. Define Parameters ---
    V_RF = 1.0  # V
    f1 = 915e6  # Hz
    Cc = 1e-12  # F
    CL = 5e-12  # F
    RL = 8e3    # Ohm
    R0 = 50.0   # Ohm
    f0 = 915e6  # Hz
    C_parasitic = 2e-15  # F

    # --- 2. Assumptions for the Model ---
    # Rectifier input resistance approximation
    R_in = RL / 2.0
    # Effective input capacitance from coupling capacitors
    C_in_rectifier = Cc / 2.0
    # Total capacitance in parallel at the rectifier input
    C_load_total = C_in_rectifier + C_parasitic
    # Power Conversion Efficiency (PCE) assumption
    PCE = 1.0

    # --- 3. Harmonic Analysis ---
    harmonics = [1, 3, 5, 7]
    V_s_amplitudes = {
        1: V_RF,
        3: V_RF * 0.9,
        5: V_RF * 0.9**2,
        7: V_RF * 0.9**3
    }
    
    V_in_sq_sum = 0.0
    v_in_magnitudes = []

    print("Calculating voltage at rectifier input for each harmonic:")
    
    for n in harmonics:
        f_n = n * f1
        w_n = 2 * math.pi * f_n
        
        # Frequency-dependent parasitic resistance
        R_p_n = R0 * (f_n / f0)**2
        
        # Load admittance (parallel R_in and C_load_total)
        Y_load_n = 1/R_in + 1j * w_n * C_load_total
        Z_load_n = 1 / Y_load_n
        
        # Source voltage for the harmonic
        V_s_n = V_s_amplitudes[n]
        
        # Voltage divider to find voltage at rectifier input
        V_in_n = V_s_n * Z_load_n / (R_p_n + Z_load_n)
        
        v_in_mag = abs(V_in_n)
        v_in_magnitudes.append(v_in_mag)
        V_in_sq_sum += v_in_mag**2
        
        print(f"  - For {n}f1 ({int(f_n/1e6)} MHz), |V_in_{n}| = {v_in_mag:.4f} V")

    # --- 4. Calculate Final DC Voltage ---
    # Based on power conservation (P_out = PCE * P_in) and our model (R_in = RL/2)
    # this simplifies to V_DC = sqrt(PCE * sum(|V_in,n|^2))
    V_DC = math.sqrt(PCE * V_in_sq_sum)
    
    # --- 5. Output the Final Equation and Result ---
    print("\nThe final DC voltage is calculated using the formula:")
    equation_str = "V_DC = sqrt("
    for i, mag in enumerate(v_in_magnitudes):
        equation_str += f"{mag:.4f}^2"
        if i < len(v_in_magnitudes) - 1:
            equation_str += " + "
    equation_str += ")"
    print(equation_str)
    
    final_eq_str = f"V_DC = sqrt({v_in_magnitudes[0]**2:.4f} + {v_in_magnitudes[1]**2:.4f} + {v_in_magnitudes[2]**2:.4f} + {v_in_magnitudes[3]**2:.4f})"
    print(final_eq_str)
    
    final_sum_str = f"V_DC = sqrt({V_in_sq_sum:.4f})"
    print(final_sum_str)
    
    print(f"\nTotal Voltage across C_L (V_DC) = {V_DC:.4f} V")
    return V_DC

if __name__ == '__main__':
    final_voltage = calculate_output_voltage()
    # The final answer format as requested.
    # print(f"<<<{final_voltage:.4f}>>>")