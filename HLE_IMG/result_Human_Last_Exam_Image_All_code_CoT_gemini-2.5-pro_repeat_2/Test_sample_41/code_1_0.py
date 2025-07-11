import math

# --- Given and Assumed Parameters ---
P_dc_inverter = 254.97  # MW, from diagram at inverter output
V_dc = 500.0            # kV, HVDC voltage
L_sys = 0.1             # H, Assumed AC system Thevenin inductance
voltage_drop_percent = 2.5 # %
harm_dist_5th = 10.0      # %
harm_dist_3rd = 5.0       # %

# --- Assumptions ---
f = 50.0                # Hz, standard system frequency
gamma_deg = 18.0          # degrees, typical inverter extinction angle

print("--- Calculation of Total Reactive Power Compensation ---")
print("\nThis calculation is based on the following key assumptions:")
print(f"1. The AC system frequency is {f} Hz.")
print(f"2. The inverter operates at a constant extinction angle (gamma) of {gamma_deg} degrees.")
print(f"3. The given line inductance L = {L_sys} H represents the Thevenin equivalent inductance of the AC system at Bus 6.")
print("4. The harmonic distortion percentages directly correspond to the required reactive power compensation as a fraction of the DC power.")
print("-" * 50)

# --- Step 1: Calculate Reactive Power Compensation for Harmonics (Q_harm) ---
print("\nStep 1: Calculate Reactive Power Compensation for Harmonics (Q_harm)")

total_harm_dist_percent = harm_dist_5th + harm_dist_3rd
q_harm = (total_harm_dist_percent / 100.0) * P_dc_inverter

print(f"The active power at the inverter is P_dc = {P_dc_inverter:.2f} MW.")
print(f"The total harmonic distortion percentage is {harm_dist_5th:.1f}% (5th) + {harm_dist_3rd:.1f}% (3rd) = {total_harm_dist_percent:.1f}%.")
print(f"Q_harm = ({total_harm_dist_percent:.1f} / 100) * {P_dc_inverter:.2f} MW")
print(f"Q_harm = {q_harm:.2f} MVar")

# --- Step 2: Calculate Reactive Power Compensation for Voltage Support (Q_volt) ---
print("\nStep 2: Calculate Reactive Power Compensation for Voltage Support (Q_volt)")

# Estimate AC voltage V_ac from V_dc
gamma_rad = math.radians(gamma_deg)
# Using the formula V_dc = 1.35 * V_ac * cos(gamma)
v_ac_kv = V_dc / (1.35 * math.cos(gamma_rad))
v_ac = v_ac_kv * 1000

print(f"\nFirst, we estimate the AC side voltage (V_ac) from the DC voltage (V_dc = {V_dc:.1f} kV).")
print(f"Using V_dc ≈ 1.35 * V_ac * cos(γ), with γ = {gamma_deg}°:")
print(f"V_ac ≈ {V_dc:.1f} kV / (1.35 * cos({gamma_deg}°)) = {v_ac_kv:.2f} kV")

# Calculate Thevenin reactance X_th
x_th = 2 * math.pi * f * L_sys

print(f"\nNext, we calculate the AC system's Thevenin reactance (X_th) at f = {f:.1f} Hz.")
print(f"X_th = 2 * π * f * L = 2 * π * {f:.1f} * {L_sys:.1f} = {x_th:.2f} Ω")

# Calculate Short Circuit Capacity S_sc
s_sc_mva = (v_ac**2) / x_th / 1e6

print(f"\nNow, we calculate the Short Circuit Capacity (S_sc) at Bus 6.")
print(f"S_sc = V_ac^2 / X_th = ({v_ac_kv:.2f} kV)^2 / {x_th:.2f} Ω = {s_sc_mva:.2f} MVA")

# Calculate Q_volt
q_volt = (voltage_drop_percent / 100.0) * s_sc_mva

print(f"\nFinally, we calculate the reactive power needed to correct the {voltage_drop_percent:.1f}% voltage drop.")
print(f"Q_volt = (Voltage Drop %) * S_sc = ({voltage_drop_percent:.1f} / 100) * {s_sc_mva:.2f} MVA")
print(f"Q_volt = {q_volt:.2f} MVar")

# --- Step 3: Calculate Total Reactive Power Compensation (Q_total) ---
print("\nStep 3: Calculate Total Reactive Power Compensation (Q_total)")

q_total = q_harm + q_volt

print("The total compensation is the sum of the harmonic and voltage support components.")
print(f"Q_total = Q_harm + Q_volt")
print(f"Q_total = {q_harm:.2f} MVar + {q_volt:.2f} MVar")
print(f"Q_total = {q_total:.2f} MVar")

print("-" * 50)
print(f"The total reactive power compensation required is {q_total:.2f} MVar.")
