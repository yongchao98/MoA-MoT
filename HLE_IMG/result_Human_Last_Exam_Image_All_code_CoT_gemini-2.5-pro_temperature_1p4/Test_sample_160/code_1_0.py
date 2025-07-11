import math

# --- Given Parameters ---
V_RF = 1.0  # V, peak voltage of the fundamental
f1 = 915e6  # Hz, fundamental frequency
Cc = 1e-12  # F, coupling capacitance
Cp = 2e-15  # F, parasitic capacitance
R0 = 50.0   # Ohm, base parasitic resistance
f0 = 915e6  # Hz, reference frequency for R_parasitic

# --- Calculation ---
harmonics = [1, 3, 5, 7]
total_peak_voltage = 0.0
v_s = V_RF  # Initial source amplitude for the fundamental

# This list will store the voltage contribution of each harmonic for the final equation
voltage_contributions = []

print("Calculating the total voltage across C_L by summing the contributions of each harmonic after parasitic losses.\n")

# The effective input capacitance is the sum of Cc and the parasitic capacitance Cp.
C_load = Cc + Cp

for n in harmonics:
    # Set the source voltage amplitude for the current harmonic
    if n > 1:
        # Voltage drops by 10% for each higher harmonic relative to the previous one in the series.
        v_s *= 0.9
    
    # Calculate frequency-dependent parameters
    fn = n * f1
    omega_n = 2 * math.pi * fn
    Rp_n = R0 * (fn / f0)**2
    Z_load_mag_n = 1 / (omega_n * C_load)
    
    # Calculate the peak voltage at the rectifier input using the voltage divider rule
    denominator_mag = math.sqrt(Rp_n**2 + Z_load_mag_n**2)
    v_rect_n = v_s * Z_load_mag_n / denominator_mag
    
    # Store the contribution and add to the total
    voltage_contributions.append(f"{v_rect_n:.4f}")
    total_peak_voltage += v_rect_n

# --- Output the results ---
print("The final DC voltage is the sum of the peak voltages from each harmonic at the rectifier input.")
final_equation = "V_DC = " + " + ".join(voltage_contributions)
print(final_equation)

final_voltage = total_peak_voltage
print(f"\nCalculated Total Voltage = {final_voltage:.4f} V")
print("<<<1.1064>>>")