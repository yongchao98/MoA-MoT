import math

# --- Given Data ---
# Active power for each load in MW
P_L1 = 1.85
P_L2 = 1.7
P_L3 = 1.75
P_L4 = 1.9
P_L5 = 2.4

# Power factor of the loads
power_factor = 0.95

# --- Step-by-step Calculation ---

# Step 1: Calculate the total active power demand (P_total)
# This is the sum of all individual loads.
P_total = P_L1 + P_L2 + P_L3 + P_L4 + P_L5

print("Step 1: Calculate the total active power (P_total).")
print(f"P_total = {P_L1} MW + {P_L2} MW + {P_L3} MW + {P_L4} MW + {P_L5} MW = {P_total:.2f} MW\n")

# Step 2: Calculate the required reactive power compensation (Q_comp)
# To mitigate the voltage drop, the ESS must supply the total reactive power demanded by the loads.
# The relationship is Q = P * tan(phi), where phi = arccos(PF).
print("Step 2: Calculate the required reactive power (Q_comp).")

# Calculate phi, the power factor angle
phi = math.acos(power_factor)

# Calculate the required reactive power compensation
Q_comp = P_total * math.tan(phi)

# Print the final equation showing all the numbers used in the calculation
print("The equation for reactive power is: Q_required = P_total * tan(arccos(PF))")
print(f"Q_required = {P_total:.2f} * tan(arccos({power_factor}))")
print(f"Q_required = {P_total:.2f} * tan({math.degrees(phi):.2f} degrees)")
print(f"Q_required = {P_total:.2f} * {math.tan(phi):.4f}")
print(f"Q_required = {Q_comp:.2f} Mvar\n")

print(f"The required reactive power compensation from the ESS is {Q_comp:.2f} Mvar.")
