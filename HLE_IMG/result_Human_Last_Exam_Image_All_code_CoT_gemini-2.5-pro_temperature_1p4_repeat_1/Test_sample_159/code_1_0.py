import math

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # --- Given Parameters ---
    V_rf_peak = 1.0  # V, peak voltage of the fundamental
    f0 = 915e6  # Hz, fundamental frequency
    R0_parasitic = 50.0  # Ohms, base parasitic resistance
    # R_load = 8000.0  # Ohms, load resistance
    # C_parasitic = 2e-15 # Farads, parasitic capacitance (ignored as its impedance is very high)
    
    # Assumption: The rectifier input impedance is matched to the characteristic impedance
    R_rect = 50.0 # Ohms

    harmonics = [1, 3, 5, 7]
    
    # --- Calculations ---
    V_peaks = {}
    V_rms = {}
    R_parasitics = {}
    Z_totals = {}
    I_rms = {}
    P_ins = {}
    P_rects = {}

    print("Step-by-step calculation per harmonic:\n")
    
    # Calculate harmonic properties
    for n in harmonics:
        if n == 1:
            V_peaks[n] = V_rf_peak
        else:
            # Voltage drops by 10% (to 90%) relative to the previous harmonic
            V_peaks[n] = V_peaks[n-2] * 0.9

        V_rms[n] = V_peaks[n] / math.sqrt(2)
        
        # Parasitic resistance R_p(f) = R0 * (f/f0)^2 = R0 * n^2
        R_parasitics[n] = R0_parasitic * (n**2)
        
        # Total series impedance
        Z_totals[n] = R_parasitics[n] + R_rect
        
        # Current for the harmonic
        I_rms[n] = V_rms[n] / Z_totals[n]
        
        # Input power from the source for the harmonic
        P_ins[n] = I_rms[n]**2 * Z_totals[n]
        
        # Power delivered to the rectifier for the harmonic
        P_rects[n] = I_rms[n]**2 * R_rect

        print(f"--- Harmonic n={n} ---")
        print(f"  V_peak          = {V_peaks[n]:.4f} V")
        print(f"  R_parasitic     = {R_parasitics[n]:.1f} Ω")
        print(f"  Z_total         = {Z_totals[n]:.1f} Ω")
        print(f"  P_in (from source) = {P_ins[n]*1000:.4f} mW")
        print(f"  P_rect (to rectifier) = {P_rects[n]*1000:.4f} mW\n")

    # Sum powers to get totals
    P_in_total = sum(P_ins.values())
    P_rect_total = sum(P_rects.values())

    # Calculate overall efficiency
    efficiency = P_rect_total / P_in_total
    
    # --- Final Output ---
    print("--- Total Power Calculation ---")
    print("The efficiency is the ratio of total power delivered to the rectifier to the total power from the source.")
    
    p_rect_sum_str = " + ".join([f"{P_rects[n]*1000:.4f}" for n in harmonics])
    p_in_sum_str = " + ".join([f"{P_ins[n]*1000:.4f}" for n in harmonics])
    
    print(f"\nEfficiency = P_rect_total / P_in_total")
    print(f"Efficiency = ({p_rect_sum_str}) mW / ({p_in_sum_str}) mW")
    print(f"Efficiency = {P_rect_total*1000:.4f} mW / {P_in_total*1000:.4f} mW")
    print(f"Efficiency = {efficiency:.4f}")
    print(f"\nOverall System Efficiency: {efficiency * 100:.2f}%")
    
    return efficiency

# Run the calculation and store the result
final_efficiency = calculate_efficiency()

# The final answer in the requested format
# Since the question asks for the efficiency value, we provide it.
# We will format it as a percentage rounded to one decimal place.
# final_efficiency is approx 0.4204, so 42.0%
final_answer = round(final_efficiency * 100, 1)
print(f"\n<<<__{final_answer}__>>>")
