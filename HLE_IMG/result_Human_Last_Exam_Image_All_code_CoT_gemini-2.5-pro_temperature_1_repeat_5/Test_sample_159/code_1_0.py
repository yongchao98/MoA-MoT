import numpy as np

# Step 1: Define constants from the problem description
V_RF_fundamental_peak = 1.0  # V
f0 = 915e6  # Hz (from ω = 2π×915MHz)
R_L = 8e3  # Ω (8kΩ)
R0_parasitic = 50.0  # Ω
C_parasitic = 2e-15  # F (2fF)
harmonic_voltage_factor = 0.9 # Voltage is 90% of the previous harmonic's voltage

# Step 2: Define harmonics and their voltages
harmonics = np.array([1, 3, 5, 7])
V_peaks = np.zeros_like(harmonics, dtype=float)
V_peaks[0] = V_RF_fundamental_peak
for i in range(1, len(harmonics)):
    V_peaks[i] = V_peaks[i-1] * harmonic_voltage_factor

# Step 3: Assume rectifier input impedance
# For a voltage doubler rectifier, a common approximation is R_in ≈ R_L / 8.
R_rect = R_L / 8.0
print(f"Key Assumption: The rectifier has a purely resistive input impedance of R_rect = R_L / 8 = {R_rect:.0f} Ω.\n")

# Step 4-6: Calculate powers for each harmonic
total_input_power = 0.0
total_rectifier_power = 0.0

print("Power calculation per harmonic:")
print("-" * 65)
print(f"{'Harmonic (n)':<15} {'V_peak (V)':<15} {'P_in (uW)':<15} {'P_out (uW)':<15}")
print("-" * 65)

for i, n in enumerate(harmonics):
    V_peak = V_peaks[i]
    f_n = n * f0
    w_n = 2 * np.pi * f_n

    # Parasitic impedance components at harmonic frequency f_n
    R_para_n = R0_parasitic * (n**2)
    XC_para_n = -1 / (w_n * C_parasitic)

    # Total input impedance at f_n
    Z_in_real = R_para_n + R_rect
    Z_in_imag = XC_para_n
    Z_in_magnitude = np.sqrt(Z_in_real**2 + Z_in_imag**2)

    # RMS current squared at f_n
    I_rms_sq = (V_peak / Z_in_magnitude)**2 / 2

    # Power calculations for harmonic f_n
    P_in_n = I_rms_sq * Z_in_real
    # Assuming ideal rectifier, P_out is the power delivered to the rectifier
    P_rect_n = I_rms_sq * R_rect

    total_input_power += P_in_n
    total_rectifier_power += P_rect_n
    
    print(f"{n:<15} {V_peak:<15.3f} {P_in_n*1e6:<15.4f} {P_rect_n*1e6:<15.4f}")

print("-" * 65)

# Step 6: Calculate total power and efficiency
P_out = total_rectifier_power
P_in = total_input_power
efficiency = (P_out / P_in) * 100 if P_in != 0 else 0

# Print final results
print("\nSummary of Total Power:")
print(f"Total Input Power (P_in) = {P_in * 1e6:.4f} µW")
print(f"Total Output Power (P_out) = {P_out * 1e6:.4f} µW")

print("\nOverall System Efficiency Calculation:")
print(f"Efficiency η = P_out / P_in")
print(f"η = {P_out * 1e6:.4f} µW / {P_in * 1e6:.4f} µW = {efficiency:.2f}%")