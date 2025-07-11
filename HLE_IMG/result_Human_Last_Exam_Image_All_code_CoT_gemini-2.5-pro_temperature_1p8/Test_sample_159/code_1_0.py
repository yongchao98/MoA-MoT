import math

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # System Parameters from the problem description
    R0 = 50.0  # Base parasitic resistance in Ohms
    R_rect = 50.0  # Assumed rectifier input resistance in Ohms
    V_fundamental_peak = 1.0  # Peak voltage of the fundamental harmonic in Volts
    
    harmonics = [1, 3, 5, 7]
    voltages = []
    current_v = V_fundamental_peak
    for _ in harmonics:
        voltages.append(current_v)
        current_v *= 0.9

    total_input_power = 0
    total_output_power = 0
    
    p_in_contributions = []
    p_out_contributions = []

    print("--- Calculating Power for Each Harmonic ---")
    for i, n in enumerate(harmonics):
        V_peak = voltages[i]
        V_rms = V_peak / math.sqrt(2)
        
        # Calculate parasitic resistance for the nth harmonic
        R_p = R0 * (n**2)
        
        # Total resistance seen by the source for this harmonic
        R_total = R_p + R_rect
        
        # Input power drawn from the source for this harmonic
        # P_in = V_rms^2 / R_total
        p_in_h = (V_rms**2) / R_total
        p_in_contributions.append(p_in_h)
        total_input_power += p_in_h
        
        # Output power delivered to the rectifier for this harmonic
        # P_out = I_rms^2 * R_rect = (V_rms / R_total)^2 * R_rect
        p_out_h = ((V_rms / R_total)**2) * R_rect
        p_out_contributions.append(p_out_h)
        total_output_power += p_out_h
        
        print(f"\nHarmonic n={n}:")
        print(f"  V_peak = {V_peak:.3f} V")
        print(f"  R_parasitic = {R0:.1f} * {n}^2 = {R_p:.1f} Ohms")
        print(f"  R_total = {R_p:.1f} + {R_rect:.1f} = {R_total:.1f} Ohms")
        print(f"  P_in({n}) = ({V_peak:.3f}/sqrt(2))^2 / {R_total:.1f} = {p_in_h:.6f} W")
        print(f"  P_out({n}) = P_in({n}) * ({R_rect:.1f} / {R_total:.1f}) = {p_out_h:.6f} W")

    # Calculate overall efficiency
    efficiency = total_output_power / total_input_power

    print("\n--- Total Power Calculation ---")
    p_in_str = " + ".join([f"{p:.6f}" for p in p_in_contributions])
    print(f"Total Input Power P_in = {p_in_str} = {total_input_power:.6f} W")
    
    p_out_str = " + ".join([f"{p:.6f}" for p in p_out_contributions])
    print(f"Total Output Power P_out = {p_out_str} = {total_output_power:.6f} W")

    print("\n--- Overall System Efficiency ---")
    print(f"Efficiency = P_out_total / P_in_total")
    print(f"Efficiency = {total_output_power:.6f} W / {total_input_power:.6f} W = {efficiency:.4f}")
    print(f"Efficiency = {efficiency * 100:.2f}%")
    
    return efficiency

# Run the calculation and store the final answer
final_efficiency = calculate_efficiency()

# Final answer format
# print(f"<<<{final_efficiency:.4f}>>>")