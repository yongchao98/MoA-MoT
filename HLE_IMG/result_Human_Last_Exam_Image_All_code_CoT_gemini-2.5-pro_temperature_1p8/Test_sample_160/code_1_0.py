import cmath
import math

def calculate_total_voltage():
    # 1. Define input parameters from the problem description
    V_RF = 1.0  # V, peak amplitude of the fundamental
    f0 = 915e6  # Hz, fundamental frequency
    Cc = 1e-12  # F, coupling capacitance
    RL = 8e3    # Ohm, load resistance
    R0 = 50     # Ohm, base parasitic resistance
    # Parasitic and load capacitances are considered part of the output filter,
    # assumed ideal in this DC model.
    # CL = 5e-12 F
    # C_parasitic = 2e-15 F

    # Harmonics to consider
    harmonics = [1, 3, 5, 7]
    
    # 2. Define harmonic properties
    V_amplitudes = {}
    V_amplitudes[1] = V_RF
    # Voltage drops by 10% for each higher harmonic relative to the previous
    for i in range(1, len(harmonics)):
        prev_h = harmonics[i-1]
        curr_h = harmonics[i]
        V_amplitudes[curr_h] = V_amplitudes[prev_h] * 0.9

    R_parasitic = {n: R0 * n**2 for n in harmonics}
    
    # 3. Model the circuit
    # The rectifier is a 2-stage voltage doubler.
    # Approximate input resistance of the rectifier
    N_stages = 2
    R_in_rect = RL / (N_stages**2)

    V_dc_components_sq = []

    print("--- Calculation Steps ---")
    print(f"Rectifier Input Resistance (R_in_rect) = R_L / N^2 = {RL:.0f} / {N_stages}^2 = {R_in_rect:.0f} Ohms\n")

    for n in harmonics:
        # Angular frequency for the harmonic
        omega_n = 2 * math.pi * f0 * n
        
        # Impedance of the two series coupling capacitors
        Z_Cc_total = 1 / (1j * omega_n * Cc / 2)
        
        # Total impedance in the divider calculation
        Z_total = R_parasitic[n] + Z_Cc_total + R_in_rect
        
        # 4. Calculate the peak voltage at the rectifier input for this harmonic
        # Using the voltage divider formula: V_out = V_in * (Z_load / Z_total)
        V_in_rect_peak = V_amplitudes[n] * abs(R_in_rect / Z_total)
        
        # The DC voltage produced by this harmonic (voltage doubler)
        V_dc_n = 2 * V_in_rect_peak
        
        # Store the square of the DC component
        V_dc_components_sq.append(V_dc_n**2)

        print(f"--- Harmonic n={n} ---")
        print(f"  Input Voltage Amplitude (V_{n}): {V_amplitudes[n]:.3f} V")
        print(f"  Parasitic Resistance (R_p({n})): {R_parasitic[n]:.1f} Ohms")
        print(f"  Coupling Capacitance Impedance (Z_Cc({n})): {Z_Cc_total.imag:.2f}j Ohms")
        print(f"  Peak AC Voltage at Rectifier Input (V_in_rect_peak({n})): {V_in_rect_peak:.4f} V")
        print(f"  DC Voltage Contribution (V_dc({n})): {V_dc_n:.4f} V")
        print("-" * 25)

    # 5. Calculate total DC voltage
    V_dc_total = math.sqrt(sum(V_dc_components_sq))

    # Print the final calculation step
    print("\n--- Final Calculation ---")
    sum_expr = " + ".join([f"{v_sq**0.5:.4f}^2" for v_sq in V_dc_components_sq])
    print(f"V_total^2 = {sum_expr}")
    sum_vals_expr = " + ".join([f"{v_sq:.4f}" for v_sq in V_dc_components_sq])
    print(f"V_total^2 = {sum_vals_expr} = {sum(V_dc_components_sq):.4f}")
    print(f"V_total = sqrt({sum(V_dc_components_sq):.4f}) = {V_dc_total:.4f} V")
    
    # Final formatted answer
    print("\nThe final calculated total voltage across C_L is:")
    print(f"{V_dc_total:.3f} V")
    return V_dc_total

if __name__ == '__main__':
    final_voltage = calculate_total_voltage()
    # The final answer in the required format
    # print(f"<<<{final_voltage:.3f}>>>") # Let's format to one decimal place as in example
    print(f"\n<<<2.7>>>")
