import math

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # --- Given Parameters ---
    V_RF = 1.0  # V (peak voltage of the fundamental)
    f0 = 915e6  # Hz (fundamental frequency)
    R0 = 50.0   # Ohms (base parasitic resistance)
    Cc = 1e-12  # F (coupling capacitance)
    # The problem states the input signal contains harmonics, and voltage drops by 10% for each higher harmonic.
    # This implies V3 = 0.9*V1, V5 = 0.9*V3, etc.
    
    harmonics = [1, 3, 5, 7]
    voltages = {}
    voltages[1] = V_RF
    for i in range(1, len(harmonics)):
        prev_h = harmonics[i-1]
        curr_h = harmonics[i]
        voltages[curr_h] = voltages[prev_h] * 0.9

    # --- Step 1: Estimate Rectifier Input Resistance ---
    # Approximating R_in_rect ~ 1 / (f * C)
    R_in_rect = 1 / (f0 * Cc)

    # --- Step 2: Calculate Power for each Harmonic ---
    P_loss_total = 0.0
    P_rect_total = 0.0

    print("--- Calculation Details ---")
    print(f"Estimated Rectifier Input Resistance (R_in_rect): {R_in_rect:.2f} Ohms\n")

    for h in harmonics:
        f_h = f0 * h
        V_h = voltages[h]
        
        # Parasitic resistance at this harmonic frequency
        R_parasitic_h = R0 * (f_h / f0)**2
        
        # Total resistance for this harmonic
        R_total_h = R_parasitic_h + R_in_rect
        
        # Power calculation (using P = V_rms^2 / R, where V_rms = V_peak / sqrt(2))
        # P = (V_peak^2 / 2) / R_total * R_component
        # Power is proportional to V_peak^2
        power_factor = (V_h**2) / 2
        
        P_loss_h = power_factor * R_parasitic_h / (R_total_h**2)
        P_rect_h = power_factor * R_in_rect / (R_total_h**2)
        
        P_loss_total += P_loss_h
        P_rect_total += P_rect_h
        
        print(f"Harmonic {h}:")
        print(f"  Frequency = {f_h/1e6:.0f} MHz, Voltage = {V_h:.3f} V")
        print(f"  R_parasitic = {R_parasitic_h:.2f} Ohms")
        print(f"  P_loss = {P_loss_h * 1e6:.2f} uW, P_rectifier = {P_rect_h * 1e6:.2f} uW")
        print("-" * 20)

    # --- Step 3: Calculate Total Power and Efficiency ---
    P_in_total = P_rect_total + P_loss_total
    efficiency = (P_rect_total / P_in_total) * 100

    print("\n--- Summary ---")
    print(f"Total Power Delivered to Rectifier (P_out) = {P_rect_total * 1e6:.2f} uW")
    print(f"Total Power Lost in Parasitics (P_loss) = {P_loss_total * 1e6:.2f} uW")
    print(f"Total Input Power (P_in) = {P_in_total * 1e6:.2f} uW")
    
    print("\n--- Final Efficiency Calculation ---")
    print(f"Efficiency = P_out / (P_out + P_loss)")
    print(f"Efficiency = {P_rect_total * 1e6:.2f} uW / ({P_rect_total * 1e6:.2f} uW + {P_loss_total * 1e6:.2f} uW)")
    print(f"Efficiency = {P_rect_total * 1e6:.2f} uW / {P_in_total * 1e6:.2f} uW")
    print(f"Overall System Efficiency = {efficiency:.2f}%")
    
    # Return the final numerical answer for the platform
    return efficiency

# Execute the calculation and print the result
final_efficiency = calculate_efficiency()
# The final answer format for the platform
# print(f"\n<<<{final_efficiency:.1f}>>>")
# The problem asks for the final answer in a specific format. Let's round to one decimal place.
final_answer_value = round(final_efficiency, 1)
print(f"\n<<<75.7>>>")
