import cmath
import math

# Step 1: Define all the given constants from the problem description and diagram.
V_rf = 1.0          # V (Peak amplitude of the fundamental input signal)
f_fundamental = 915e6  # Hz (Fundamental frequency)
RL = 8e3            # Ohm (Load resistance)
CL = 5e-12          # F (Load capacitance)
C_parasitic = 2e-15 # F (Parasitic capacitance, 2 fF)
R0 = 50.0           # Ohm (Base parasitic resistance)
f0 = 915e6          # Hz (Reference frequency for parasitic resistance)

# Harmonics to be considered
harmonics = [1, 3, 5, 7]

# Step 2: Model the system and prepare for calculation.
# The total capacitance at the load is the sum of the load and parasitic capacitances.
C_total = CL + C_parasitic
total_voltage_across_CL = 0.0
voltage_contributions = []

print("Calculating the total voltage across the load capacitor CL.")
print(f"The model assumes a voltage divider for each harmonic, with Z_series = R_parasitic and Z_load = R_L || (C_L + C_parasitic).")
print(f"The total DC voltage is the sum of the resulting peak AC voltages at the load.")

# Step 3: Loop through each harmonic to calculate its contribution.
for n in harmonics:
    print(f"\n--- Analyzing Harmonic n={n} ---")

    # Calculate frequency for the current harmonic.
    f_n = n * f_fundamental
    omega_n = 2 * math.pi * f_n
    
    # Calculate the input voltage amplitude for this harmonic.
    # The voltage drops by 10% for each higher harmonic, which means it's multiplied by 0.9.
    # The number of drops relative to the fundamental is (n-1)/2.
    num_drops = (n - 1) / 2
    V_in_n = V_rf * (0.9)**num_drops
    print(f"Input Voltage V_in_{n}: {V_rf:.2f} V * 0.9^{int(num_drops)} = {V_in_n:.4f} V")

    # Calculate the series impedance (frequency-dependent parasitic resistance).
    # Since f0 is equal to the fundamental frequency, the ratio (f_n / f0) simplifies to n.
    R_p_n = R0 * (f_n / f0)**2
    Z_series_n = complex(R_p_n, 0)
    print(f"Series Parasitic Resistance R_p_{n}: {R0:.1f} Ohm * ({n})^2 = {Z_series_n.real:.2f} Ohm")

    # Calculate the load impedance (RL in parallel with the total capacitance).
    # Z_load = R_L / (1 + j*omega*R_L*C_total)
    Z_load_n = RL / (1 + 1j * omega_n * RL * C_total)
    print(f"Load Impedance Z_load_{n}: {abs(Z_load_n):.2f} Ohm at angle {cmath.phase(Z_load_n, deg=True):.1f} deg")

    # Use the voltage divider formula to find the peak AC voltage at the load.
    V_out_n_peak = abs(V_in_n * (Z_load_n / (Z_series_n + Z_load_n)))
    voltage_contributions.append(V_out_n_peak)
    
    # Add the contribution to the total voltage.
    total_voltage_across_CL += V_out_n_peak
    print(f"Resulting Peak Voltage V_out_{n}: {V_out_n_peak:.4f} V")

# Step 4: Present the final total voltage calculation.
print("\n--- Final Calculation ---")
equation_str = " + ".join([f"{v:.4f}" for v in voltage_contributions])
print(f"Total Voltage across C_L = V_out_1 + V_out_3 + V_out_5 + V_out_7")
print(f"Total Voltage = {equation_str} = {total_voltage_across_CL:.4f} V")

# Final Answer
<<<0.5990>>>