import cmath
import math

def calculate_system_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and
    parasitic losses based on a simplified linear AC model.
    """
    # --- Givens from the problem description and table ---
    f_fund = 915e6  # Fundamental frequency in Hz (from ω = 2π*915MHz)
    V_fund = 1.0    # Fundamental peak voltage in V
    R_L = 8e3       # Load resistance in Ohm
    C_L = 5e-12     # Load capacitance in Farad
    R0 = 50.0       # Base parasitic resistance in Ohm
    f0 = 915e6      # Base frequency for parasitic resistance in Hz
    C_parasitic = 2e-15 # Parasitic capacitance in Farad

    # The model assumes the parasitic capacitance is in parallel with the main load.
    C_total_parallel = C_L + C_parasitic

    # --- Harmonics Information ---
    harmonics = [1, 3, 5, 7]
    voltages = {1: V_fund}
    # Voltage drops by 10% for each higher harmonic relative to the previous
    for i in range(1, len(harmonics)):
        prev_n = harmonics[i-1]
        curr_n = harmonics[i]
        voltages[curr_n] = voltages[prev_n] * 0.9

    total_input_power = 0
    total_output_power = 0

    print("Analyzing the system for each harmonic based on a simplified AC model:")
    print("="*60)

    for n in harmonics:
        f_n = n * f_fund
        w_n = 2 * math.pi * f_n
        V_in_n = voltages[n]

        # 1. Calculate parasitic resistance for this harmonic
        R_p_n = R0 * (f_n / f0)**2
        
        # 2. Calculate the combined load impedance Z_load = R_L || (C_L + C_parasitic)
        # Admittance of the parallel components: Y_load = 1/R_L + j*w*C_total
        Y_load_n = 1/R_L + 1j * w_n * C_total_parallel
        Z_load_n = 1 / Y_load_n
        
        # 3. Calculate total impedance seen by the source: Z_total = R_p + Z_load
        Z_total_n = R_p_n + Z_load_n
        
        # 4. Calculate input power for this harmonic
        # P_in = 0.5 * V_peak^2 * Re(Z_total) / |Z_total|^2
        P_in_n = 0.5 * (V_in_n**2) * Z_total_n.real / abs(Z_total_n)**2
        total_input_power += P_in_n
        
        # 5. Calculate output power (power dissipated in R_L) for this harmonic
        # Use voltage divider to find voltage across the load: V_load = V_in * Z_load / Z_total
        # Power in R_L: P_out = 0.5 * |V_load|^2 / R_L
        V_load_n = V_in_n * Z_load_n / Z_total_n
        P_out_n = 0.5 * abs(V_load_n)**2 / R_L
        total_output_power += P_out_n
        
        print(f"Harmonic n={n}:")
        print(f"  - Frequency: {f_n/1e6:.0f} MHz")
        print(f"  - Input Voltage: {V_in_n:.3f} V")
        print(f"  - Parasitic Resistance: {R_p_n:.1f} Ohm")
        print(f"  - Combined Load Impedance: {Z_load_n.real:.3f} + {Z_load_n.imag:.3f}j Ohm")
        print(f"  - Total System Impedance: {Z_total_n.real:.3f} + {Z_total_n.imag:.3f}j Ohm")
        print(f"  - Input Power: {P_in_n * 1000:.4f} mW")
        print(f"  - Output Power (in R_L): {P_out_n * 1e6:.4f} uW")
        print("-" * 25)

    # --- Final Calculation ---
    # Avoid division by zero if input power is zero
    efficiency = total_output_power / total_input_power if total_input_power > 0 else 0

    print("="*60)
    print("Summary:")
    print(f"Total Input Power (sum of all harmonics): {total_input_power * 1000:.4f} mW")
    print(f"Total Output Power (sum of power in R_L): {total_output_power * 1000:.4f} mW")
    print("\nOverall System Efficiency Calculation:")
    print(f"Efficiency η = Total Output Power / Total Input Power")
    # Using .6g format to show numbers concisely
    print(f"η = {total_output_power:.6g} W / {total_input_power:.6g} W")
    print(f"η = {efficiency:.6f}")
    print(f"η = {efficiency * 100:.4f} %")
    
    # Returning the final numerical answer for the <<<>>> format
    return efficiency

# Run the calculation and store the result
final_efficiency = calculate_system_efficiency()

# The final answer is printed in the required format
# print(f"\n<<<0.002540>>>")
# Using the calculated value for better precision
print(f"\n<<<{final_efficiency:.6f}>>>")