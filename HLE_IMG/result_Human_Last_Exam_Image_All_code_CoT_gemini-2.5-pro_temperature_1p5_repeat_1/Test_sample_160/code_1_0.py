import math
import cmath

# --- Define constants from the problem description and circuit diagram ---
V_RF = 1.0  # V, peak voltage of the fundamental signal
f1 = 915e6  # Hz, fundamental frequency
R_L = 8000  # Ohms, load resistance
C_L = 5e-12  # Farads, load capacitance

# Parasitic components
C_parasitic = 2e-15  # Farads (2fF)
R0 = 50  # Ohms, base parasitic resistance
f0 = 915e6  # Hz, reference frequency for parasitic resistance

# --- Calculations ---

# 1. Calculate the total capacitance in parallel with the load resistor
C_total = C_L + C_parasitic

# 2. Define the harmonics to be considered
harmonics = [1, 3, 5, 7]
voltages = {}
powers = {}
total_power = 0.0

# 3. Calculate the peak voltage for each harmonic
# The voltage drops by 10% for each higher harmonic relative to the previous one
v_peak_current = V_RF
voltages[1] = v_peak_current
for i in range(len(harmonics) - 1):
    v_peak_current *= 0.90
    voltages[harmonics[i+1]] = v_peak_current

# 4. Loop through each harmonic to calculate its power contribution
print("Calculating power contribution from each harmonic:")
print("-" * 60)

for n in harmonics:
    # a. Calculate the frequency of the harmonic
    f_n = n * f1
    omega_n = 2 * math.pi * f_n

    # b. Get the input voltage for the harmonic
    v_in_n = voltages[n]

    # c. Calculate the frequency-dependent parasitic resistance
    R_p_n = R0 * (f_n / f0)**2

    # d. Calculate the complex impedance of the parallel RC load
    Z_L_n = 1 / (1/R_L + 1j * omega_n * C_total)

    # e. Calculate the total complex impedance of the series-parallel circuit
    Z_total_n = R_p_n + Z_L_n

    # f. Calculate the average power delivered to the load impedance
    # P_avg = (1/2) * |V_peak / Z_total|^2 * Re(Z_load)
    power_n = 0.5 * (v_in_n / abs(Z_total_n))**2 * Z_L_n.real
    powers[n] = power_n
    total_power += power_n

    print(f"For harmonic n={n}:")
    print(f"  f_{n:<2} = {f_n/1e6:8.2f} MHz")
    print(f"  V_{n:<2} = {v_in_n:.3f} V")
    print(f"  R_parasitic({n}) = {R_p_n:7.2f} Ohms")
    print(f"  Z_L({n})         = {Z_L_n.real:7.4f} + {Z_L_n.imag:7.4f}j Ohms")
    print(f"  Z_total({n})     = {Z_total_n.real:7.4f} + {Z_total_n.imag:7.4f}j Ohms")
    print(f"  P_L({n})         = {power_n * 1000:.4f} mW\n")

# 5. Calculate the total power delivered to the load
print("-" * 60)
print("Calculating total power and final voltage:")
p_strings = [f"{powers[n]*1000:.4f}" for n in harmonics]
print(f"P_total = P_L(1) + P_L(3) + P_L(5) + P_L(7)")
print(f"P_total = {' + '.join(p_strings)} mW")
print(f"P_total = {total_power * 1000:.4f} mW")
print("")

# 6. Calculate the final output DC voltage across the load
# Assuming 100% rectifier efficiency, P_DC = P_total_AC
# P_DC = V_DC^2 / R_L  =>  V_DC = sqrt(P_total_AC * R_L)
v_dc = math.sqrt(total_power * R_L)

print(f"Final Voltage Calculation (assuming 100% rectification efficiency):")
print(f"V_DC = sqrt(P_total * R_L)")
print(f"V_DC = sqrt({total_power:.6f} W * {R_L:.0f} Ω)")
print(f"V_DC = sqrt({total_power * R_L:.4f})")
print(f"V_DC = {v_dc:.4f} V")