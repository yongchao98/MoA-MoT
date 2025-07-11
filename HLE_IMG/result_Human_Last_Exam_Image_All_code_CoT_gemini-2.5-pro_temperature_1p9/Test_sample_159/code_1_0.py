import math

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # 1. Define constants and assumptions from the problem description
    V_fundamental_peak = 1.0  # V, peak voltage of the fundamental frequency
    f_fundamental = 915e6  # Hz
    R_L = 8e3  # Ω
    
    # Parasitic component parameters
    R0_parasitic = 50.0  # Ω
    f0_parasitic = 915e6  # Hz
    
    # Harmonic distortion parameter
    harmonic_drop_factor = 0.9
    
    # Assumption: Effective input resistance of the voltage doubler rectifier
    R_in_eff = R_L / 2
    
    harmonics = [1, 3, 5, 7]
    V_peak = V_fundamental_peak
    total_source_power = 0
    total_rectifier_power = 0

    print("--- System Parameters and Assumptions ---")
    print(f"Load Resistance (R_L): {R_L/1000} kΩ")
    print(f"Assumed Rectifier Effective Input Resistance (R_in_eff = R_L/2): {R_in_eff/1000} kΩ")
    print(f"Parasitic Resistance Model: R_parasitic(f) = {R0_parasitic} * (f / {f0_parasitic/1e6} MHz)^2 Ω\n")

    print("--- Analysis per Harmonic ---")
    for n in harmonics:
        if n > 1:
            # Voltage drops by 10% for each higher harmonic
            V_peak *= harmonic_drop_factor

        f_n = n * f_fundamental
        R_parasitic = R0_parasitic * (f_n / f0_parasitic)**2
        
        # Total resistance in the circuit for this harmonic
        R_total = R_parasitic + R_in_eff
        
        # Total power drawn from the source for this harmonic
        # Formula: P = V_rms^2 / R = (V_peak/√2)^2 / R = V_peak^2 / (2 * R)
        power_source_n = (V_peak**2) / (2 * R_total)
        
        # Power delivered to the rectifier's effective resistance
        # Formula: P = I_rms^2 * R_in_eff, where I_rms = V_rms / R_total
        # P = (V_peak^2 * R_in_eff) / (2 * R_total^2)
        power_rectifier_n = (V_peak**2 * R_in_eff) / (2 * R_total**2)
        
        total_source_power += power_source_n
        total_rectifier_power += power_rectifier_n
        
        print(f"Harmonic {n}: V_peak = {V_peak:.3f} V, R_parasitic = {R_parasitic:.1f} Ω, "
              f"P_source = {power_source_n*1e6:.2f} µW, P_delivered = {power_rectifier_n*1e6:.2f} µW")

    # Calculate overall efficiency
    efficiency = total_rectifier_power / total_source_power
    
    print("\n--- Total Power Calculation ---")
    print(f"Total Power Supplied by Source: {total_source_power * 1e6:.2f} µW")
    print(f"Total Power Delivered to Rectifier: {total_rectifier_power * 1e6:.2f} µW")
    
    print("\n--- Final Efficiency Calculation ---")
    print(f"Efficiency η = (Total Power Delivered to Rectifier) / (Total Power Supplied by Source)")
    print(f"η = {total_rectifier_power:.6f} W / {total_source_power:.6f} W")
    print(f"η = {efficiency:.4f}")
    print(f"The overall system efficiency is {efficiency * 100:.2f}%")
    
    return efficiency

# Execute the calculation
final_efficiency = calculate_efficiency()
# The final numerical answer for the platform
# Return the efficiency as a decimal value, not percentage.
print(f"<<<{final_efficiency:.4f}>>>")