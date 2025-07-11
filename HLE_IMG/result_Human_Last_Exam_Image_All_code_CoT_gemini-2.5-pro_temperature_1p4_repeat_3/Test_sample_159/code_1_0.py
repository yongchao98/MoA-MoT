import math

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # --- Given Parameters ---
    V_fund_peak = 1.0  # V, peak voltage of the fundamental frequency
    f0 = 915e6         # Hz, fundamental frequency
    R0 = 50.0          # Ohms, base parasitic resistance
    RL = 8000.0        # Ohms, load resistance
    # The problem specifies 3rd, 5th, and 7th harmonics. n=1 is the fundamental.
    harmonics = [1, 3, 5, 7]
    voltage_drop_factor = 0.9 # Voltage is 90% of the previous harmonic's voltage

    # --- Step 1: Calculate Output Power (P_out) ---
    print("--- Calculating Output Power (P_out) ---")

    # Calculate peak voltages for all specified harmonics
    V_peaks = {}
    current_V_peak = V_fund_peak
    V_peaks[1] = current_V_peak
    print(f"V1 (peak) = {V_peaks[1]:.3f} V")
    for i in range(1, len(harmonics)):
        # Voltage drops by 10% relative to the PREVIOUS harmonic
        current_V_peak *= voltage_drop_factor
        V_peaks[harmonics[i]] = current_V_peak
        print(f"V{harmonics[i]} (peak) = {V_peaks[harmonics[i]]:.3f} V")

    # The total peak voltage is the sum, as odd harmonics' peaks align
    V_dc = sum(V_peaks.values())
    print(f"\nTotal peak voltage (V_peak = V_DC) = " + " + ".join([f"{v:.3f}" for v in V_peaks.values()]) + f" = {V_dc:.3f} V")

    # Calculate DC output power
    p_out = (V_dc**2) / RL
    print(f"P_out = V_DC^2 / R_L = ({V_dc:.3f} V)^2 / {RL:.0f} Ω = {p_out:.6f} W")

    # --- Step 2: Calculate Parasitic Power Loss (P_loss) ---
    print("\n--- Calculating Parasitic Power Loss (P_loss) ---")

    # Calculate parasitic resistance and power loss for each harmonic
    p_loss_total = 0
    p_losses = {}
    r_parasitics = {}

    for n in harmonics:
        # R_parasitic(f) = R0 * (f/f0)^2. Since f = n*f0, R_p,n = R0 * n^2
        r_p_n = R0 * (n**2)
        r_parasitics[n] = r_p_n
        
        # P_loss,n = V_peak,n^2 / (2 * R_p,n)
        v_n = V_peaks[n]
        p_loss_n = (v_n**2) / (2 * r_p_n)
        p_losses[n] = p_loss_n
        p_loss_total += p_loss_n
        print(f"For n={n}: R_parasitic = {r_p_n:.0f} Ω, P_loss,{n} = ({v_n:.3f} V)^2 / (2 * {r_p_n:.0f} Ω) = {p_loss_n:.6f} W")

    print("\nTotal Parasitic Power Loss (P_loss) = " + " + ".join([f"{p:.6f}" for p in p_losses.values()]) + f" = {p_loss_total:.6f} W")
    
    # --- Step 3: Calculate Total Input Power (P_in) ---
    p_in = p_out + p_loss_total
    print(f"\n--- Calculating Total Input Power (P_in) ---")
    print(f"P_in = P_out + P_loss = {p_out:.6f} W + {p_loss_total:.6f} W = {p_in:.6f} W")

    # --- Step 4: Calculate Overall System Efficiency ---
    efficiency = (p_out / p_in)
    efficiency_percent = efficiency * 100
    print(f"\n--- Calculating Final Efficiency (η) ---")
    print(f"Efficiency (η) = P_out / P_in = {p_out:.6f} W / {p_in:.6f} W = {efficiency:.4f}")
    print(f"\nOverall system efficiency: {efficiency_percent:.2f}%")

    return efficiency_percent

if __name__ == '__main__':
    final_efficiency = calculate_efficiency()
    # The final answer is wrapped according to the instruction format.
    print(f"\n<<<System Efficiency: {final_efficiency:.2f} %>>>")

# Execute the function to get the final answer in the requested format for the user
# calculate_efficiency() # This line is commented out as the final answer will be wrapped outside the code block.
