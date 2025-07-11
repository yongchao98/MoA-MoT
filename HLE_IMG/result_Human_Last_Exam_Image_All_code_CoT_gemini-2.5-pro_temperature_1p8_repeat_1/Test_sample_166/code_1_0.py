import math

# Step 1: Define constants and given values
loads_active_power = {
    'L1': 1.85,  # MW
    'L2': 1.7,   # MW
    'L3': 1.75,  # MW
    'L4': 1.9,   # MW
    'L5': 2.4,   # MW
}
power_factor = 0.95
line_length = 10  # km
resistance_per_km = 0.03  # Ohm/km
voltage_drop_percent = 0.05
nominal_voltage_kv = 22  # kV

# Assumption: Assume a standard X/R ratio of 3 as reactance is not given.
X_R_ratio = 3

# Step 2: Calculate total active and reactive power of the loads
p_total_mw = sum(loads_active_power.values())
angle = math.acos(power_factor)
q_total_mvar = p_total_mw * math.tan(angle)

print(f"1. Calculating Total System Load:")
print(f"   - Total Active Power (P_total): {p_total_mw:.3f} MW")
print(f"   - Total Reactive Power (Q_load): {q_total_mvar:.3f} Mvar\n")

# Step 3: Calculate line impedance
r_line = resistance_per_km * line_length
x_line = r_line * X_R_ratio

print(f"2. Calculating Line Impedance:")
print(f"   - Line Resistance (R): {r_line:.3f} Ohms")
print(f"   - Assumed X/R Ratio: {X_R_ratio}")
print(f"   - Line Reactance (X): {x_line:.3f} Ohms\n")

# Step 4: Define voltage conditions
v_sending_kv = nominal_voltage_kv
v_receiving_kv = nominal_voltage_kv * (1 - voltage_drop_percent)
v_drop_kv = v_sending_kv - v_receiving_kv

print(f"3. Calculating Voltage Conditions:")
print(f"   - Sending End Voltage (V_s) at Bus 1: {v_sending_kv:.2f} kV")
print(f"   - Target Receiving End Voltage (V_r) at Bus 4 (after 5% drop): {v_receiving_kv:.2f} kV")
print(f"   - Absolute Voltage Drop (ΔV): {v_drop_kv:.2f} kV\n")

# Step 5 & 6: Use the voltage drop formula to find the net reactive power flow (Q_net)
# ΔV ≈ (R * P + X * Q_net) / V_r  (where P, Q are in MW/Mvar and V in kV)
# Rearranging for Q_net: Q_net = (ΔV * V_r - R * P) / X
q_net_mvar = (v_drop_kv * v_receiving_kv - r_line * p_total_mw) / x_line

print(f"4. Calculating Net Reactive Power Flow (Q_net):")
print(f"   - Using formula: Q_net = (ΔV * V_r - R * P_total) / X")
print(f"   - Q_net = ({v_drop_kv:.2f} kV * {v_receiving_kv:.2f} kV - {r_line:.3f} Ω * {p_total_mw:.3f} MW) / {x_line:.3f} Ω")
print(f"   - Q_net = ({v_drop_kv * v_receiving_kv:.3f} - {r_line * p_total_mw:.3f}) / {x_line:.3f}")
print(f"   - Q_net = {q_net_mvar:.3f} Mvar\n")

# Step 7: Calculate the required ESS compensation
# Q_net = Q_load - Q_ESS  (A positive Q_ESS means injection)
# Rearranging for Q_ESS: Q_ESS = Q_load - Q_net
q_ess_mvar = q_total_mvar - q_net_mvar

# The question asks for the compensation *from* the ESS, which implies injection.
# A negative result from our calculation means injection. We report the magnitude.
q_ess_compensation_magnitude = abs(q_ess_mvar)

print(f"5. Calculating Required Reactive Power from ESS (Q_ESS):")
print(f"   - The relationship is: Q_net = Q_load - Q_ESS")
print(f"   - Rearranging for Q_ESS: Q_ESS = Q_load - Q_net")
print(f"   - Final Equation: Q_ESS = {q_total_mvar:.3f} Mvar - {q_net_mvar:.3f} Mvar")
print(f"   - Q_ESS = {q_ess_mvar:.3f} Mvar\n")
print(f"The negative sign indicates that the ESS must inject reactive power.")
print(f"The required reactive power compensation from the ESS is {q_ess_compensation_magnitude:.2f} Mvar.")
print("Final numerical equation for compensation (in Mvar):")
print(f"{q_ess_compensation_magnitude:.2f} = |{q_total_mvar:.3f} - (({v_drop_kv:.2f} * {v_receiving_kv:.2f}) - ({r_line:.3f} * {p_total_mw:.3f})) / {x_line:.3f}|")

# Final Answer
print(f"\n<<<The required reactive power compensation from the ESS is {q_ess_compensation_magnitude:.2f} Mvar>>>")
