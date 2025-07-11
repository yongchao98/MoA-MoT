import math

# Step 1: Define the given data from the power distribution network diagram.
# Real power of each load in MW.
P_L1 = 1.85
P_L2 = 1.7
P_L3 = 1.75
P_L4 = 1.9
P_L5 = 2.4
# Power factor for all loads.
power_factor = 0.95

# Step 2: Calculate the total real power (P_total) by summing all loads.
P_total = P_L1 + P_L2 + P_L3 + P_L4 + P_L5

# Step 3: Calculate the required reactive power compensation (Q_ESS).
# This is equal to the total reactive power consumed by the loads (Q_total).
# The formula relating real power (P), reactive power (Q), and power factor (PF) is:
# Q = P * tan(phi), where phi = arccos(PF).
phi = math.acos(power_factor)
tan_phi = math.tan(phi)
Q_ESS = P_total * tan_phi

# Step 4: Print the detailed calculation steps and the final answer.
print("To determine the required reactive power compensation, we first sum the real power of all loads and then use the power factor to find the total reactive power demand.")
print("\n--- Calculation Steps ---")

# Print the calculation for total real power.
print("\n1. Total Real Power (P_total):")
print(f"P_total = P_L1 + P_L2 + P_L3 + P_L4 + P_L5")
print(f"P_total = {P_L1} MW + {P_L2} MW + {P_L3} MW + {P_L4} MW + {P_L5} MW")
print(f"P_total = {P_total:.2f} MW")

# Print the calculation for reactive power.
print("\n2. Required Reactive Power Compensation (Q_ESS):")
print("Q_ESS = P_total * tan(arccos(PF))")
print(f"Q_ESS = {P_total:.2f} MW * tan(arccos({power_factor}))")
print(f"Q_ESS = {P_total:.2f} MW * {tan_phi:.4f}")
print(f"Q_ESS = {Q_ESS:.2f} MVAR")

print("\n-------------------------")
print(f"The final required reactive power compensation from the ESS is {Q_ESS:.2f} MVAR.")