import math

# Step 1 & 2: Define given values and calculate the fundamental reactive power compensation.
# Active power transmitted by the inverter (P_inverter) in MW.
P_inverter = 254.97
# Net reactive power exported to the grid from the inverter station in MVar.
Q_export = 14.58
# Heuristic factor for inverter reactive power consumption (Q ≈ 0.5 * P).
inverter_q_factor = 0.5

# Calculate the reactive power consumed by the inverter.
Q_inverter = inverter_q_factor * P_inverter

# The fundamental compensation must supply the inverter's need plus the exported amount.
Q_comp_fund = Q_inverter + Q_export

# Step 3: Calculate the reactive power compensation needed for harmonics.
# Harmonic distortion percentages.
h5_distortion_percent = 10
h3_distortion_percent = 5
total_harmonic_distortion_percent = h5_distortion_percent + h3_distortion_percent

# Calculate harmonic reactive power as a percentage of active power.
Q_harmonics = (total_harmonic_distortion_percent / 100) * P_inverter

# Step 4: Calculate the total required reactive power compensation.
# Total compensation is the sum of fundamental and harmonic requirements.
total_reactive_compensation = Q_comp_fund + Q_harmonics

# Step 5: Print the final results in a clear equation format.
print("Calculation of Total Reactive Power Compensation:")
print(f"1. Fundamental compensation requirement (supplying inverter and grid export):")
print(f"   Q_inverter_consumption = {inverter_q_factor} * {P_inverter:.2f} MW = {Q_inverter:.2f} MVar")
print(f"   Q_comp_fund = Q_inverter_consumption + Q_export = {Q_inverter:.2f} MVar + {Q_export:.2f} MVar = {Q_comp_fund:.2f} MVar")
print("\n2. Harmonic reactive power requirement:")
print(f"   Q_harmonics = ({h5_distortion_percent}% + {h3_distortion_percent}%) * {P_inverter:.2f} MW = {Q_harmonics:.2f} MVar")
print("\n3. Total Reactive Power Compensation Required:")
print(f"   Total Compensation = Fundamental Requirement + Harmonic Requirement")
print(f"   Total Compensation = {Q_comp_fund:.2f} MVar + {Q_harmonics:.2f} MVar = {total_reactive_compensation:.2f} MVar")

# Final Answer
# print(f"<<<{total_reactive_compensation:.2f}>>>")