import math

# Overall plan:
# 1. Define the system parameters based on the problem description. This includes the source voltage, harmonic content,
#    parasitic elements (resistance and capacitance), and load characteristics.
# 2. Model the input voltage signal, including the fundamental frequency and the specified third, fifth, and seventh harmonics.
#    The amplitude of each harmonic is calculated based on the "voltage drops by 10% relative to the previous harmonic" rule.
# 3. Model the parasitic resistance as a frequency-dependent internal resistance of the source, R_s(f).
# 4. We will make a key simplifying assumption that the input impedance of the rectifier circuit is purely resistive
#    and equal to 50 Ohms. This assumption is based on the provided R0 = 50 Ohm value, which often implies
#    a system matched to 50 Ohms at the fundamental frequency for maximum power transfer.
# 5. We will find that the parasitic capacitance has a very high impedance compared to the rectifier load and can be neglected.
# 6. Calculate the total power delivered by the source (P_input), which is the sum of the power delivered at the fundamental
#    and each harmonic frequency.
# 7. Calculate the total power delivered to the rectifier (P_delivered), which is the sum of power dissipated in the 50 Ohm
#    rectifier load at each frequency. This represents the useful power transferred to the next stage.
# 8. The overall system efficiency is then calculated as the ratio of the total delivered power to the total input power
#    (Efficiency = P_delivered / P_input). This represents the efficiency of the power transfer from the source to the rectifier,
#    accounting for losses in the frequency-dependent parasitic resistance.

# Step 1: Define system parameters
V_rf_peak = 1.0  # V, peak voltage of the fundamental
f0 = 915e6  # Hz, fundamental frequency
R0 = 50.0  # Ohms, base parasitic resistance
C_parasitic = 2e-15  # F, parasitic capacitance
# Step 4: Assume rectifier input impedance is 50 Ohm resistive
R_rect = 50.0  # Ohms

# Step 2: Model the input signal with harmonics
harmonics = [1, 3, 5, 7]
V_peaks = {}
current_v_peak = V_rf_peak
V_peaks[1] = current_v_peak
for i in range(1, len(harmonics)):
    # Voltage drops by 10% -> new voltage is 90% of the previous harmonic's voltage
    current_v_peak *= 0.9
    V_peaks[harmonics[i]] = current_v_peak

total_input_power = 0
total_delivered_power = 0

print("Calculating efficiency based on power transfer to the rectifier.")
print("-" * 50)

# Steps 5 & 7: Loop through each harmonic to calculate power
for n in harmonics:
    f_n = n * f0
    w_n = 2 * math.pi * f_n
    
    # Check impedance of parasitic capacitance
    Z_c_parasitic_mag = 1 / (w_n * C_parasitic)
    
    # Step 3: Calculate frequency-dependent parasitic resistance
    R_parasitic_n = R0 * (f_n / f0)**2
    
    # RMS voltage squared of the source for this harmonic
    V_rms_sq_n = (V_peaks[n]**2) / 2
    
    # Total resistance for this harmonic in the simplified model
    # Model: Source -> R_parasitic_n -> R_rect
    R_total_n = R_parasitic_n + R_rect
    
    # Total power drawn from the source at this harmonic
    P_input_n = V_rms_sq_n / R_total_n
    
    # Power delivered to the rectifier load at this harmonic
    # This can be calculated using power divider rule: P_del = P_in * (R_rect / R_total)
    P_delivered_n = P_input_n * (R_rect / R_total_n)
    
    total_input_power += P_input_n
    total_delivered_power += P_delivered_n
    
    print(f"Harmonic {n}:")
    print(f"  Frequency = {f_n/1e6:.0f} MHz")
    print(f"  V_peak = {V_peaks[n]:.3f} V")
    print(f"  Parasitic Resistance (R_s) = {R_parasitic_n:.1f} Ohms")
    print(f"  Parasitic Capacitance Impedance |Z_c| = {Z_c_parasitic_mag/1e3:.1f} kOhms (Negligible)")
    print(f"  Input Power = {P_input_n * 1000:.4f} mW")
    print(f"  Delivered Power = {P_delivered_n * 1000:.4f} mW")
    print("-" * 20)

# Step 8: Calculate overall efficiency
efficiency = total_delivered_power / total_input_power

print("Final Calculation:")
print(f"Total Power Delivered to Rectifier (P_delivered) = {total_delivered_power * 1000:.4f} mW")
print(f"Total Power Drawn from Source (P_input) = {total_input_power * 1000:.4f} mW")
print("\nOverall System Efficiency (η) = P_delivered / P_input")
print(f"η = {total_delivered_power * 1000:.4f} mW / {total_input_power * 1000:.4f} mW")
print(f"η = {efficiency:.4f}")
print(f"η = {efficiency * 100:.2f}%")