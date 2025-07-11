import math

# --- Given Data ---
# Loads (in MW)
P_L1 = 1.85
P_L2 = 1.7
P_L3 = 1.75
P_L4 = 1.9
P_L5 = 2.4
# Power Factor
power_factor = 0.95

# --- Step 1: Calculate Total Real Power ---
# Summing up all the individual loads
total_real_power = P_L1 + P_L2 + P_L3 + P_L4 + P_L5

print("Step 1: Calculate the total real power demand (P_total) from all loads.")
print(f"P_total = {P_L1} MW + {P_L2} MW + {P_L3} MW + {P_L4} MW + {P_L5} MW")
print(f"P_total = {total_real_power:.2f} MW\n")

# --- Step 2: Calculate Total Reactive Power ---
# The power factor angle, phi = arccos(PF)
phi = math.acos(power_factor)
# Total reactive power, Q = P * tan(phi)
total_reactive_power = total_real_power * math.tan(phi)

print("Step 2: Calculate the total reactive power demand (Q_total) using the power factor.")
print("The formula is: Q_total = P_total * tan(arccos(Power Factor))")
print("The final equation with all the numbers is:")
# Showing the full equation with numbers
print(f"Q_total = ({P_L1} + {P_L2} + {P_L3} + {P_L4} + {P_L5}) * tan(arccos({power_factor}))")
print(f"Q_total = {total_real_power:.2f} * {math.tan(phi):.4f}")
print(f"Q_total = {total_reactive_power:.4f} MVAR\n")

# --- Step 3: Determine ESS Reactive Power Compensation ---
# The required compensation from the ESS is equal to the total reactive power demand.
Q_ESS = total_reactive_power

print("Step 3: Determine the required reactive power compensation from the ESS (Q_ESS).")
print("The ESS must supply the reactive power demanded by the loads to support the system voltage.")
print(f"Therefore, the required reactive power from the ESS is equal to the total reactive power demand.")
print(f"Q_ESS = {Q_ESS:.2f} MVAR")
