import math

# Step 1: Define the real power of each load in MW from the diagram.
L1 = 1.85
L2 = 1.7
L3 = 1.75
L4 = 1.9
L5 = 2.4

# The power factor is given for all loads.
power_factor = 0.95

# Step 2: Calculate the total real power demand by summing all loads.
total_real_power = L1 + L2 + L3 + L4 + L5

# Step 3: Calculate the total reactive power demand.
# This is the amount of reactive power the ESS needs to supply for compensation.
# The formula is Q = P * tan(acos(PF)).

# Calculate the power factor angle (phi).
phi = math.acos(power_factor)

# Calculate the total reactive power.
total_reactive_power = total_real_power * math.tan(phi)

# Step 4: Print the calculation steps and the final result.
print("Plan:")
print("1. Calculate the total real power (P_total) by summing all loads.")
print("2. Calculate the total reactive power (Q_total) using the formula: Q = P_total * tan(acos(PF)).")
print("3. The result, Q_total, is the required reactive power compensation from the ESS.\n")

print("Calculation:")
print(f"Total Real Power (P_total) = {L1} + {L2} + {L3} + {L4} + {L5} = {total_real_power:.2f} MW")
print(f"Power Factor (PF) = {power_factor}\n")

print("Final Equation:")
print(f"Required Reactive Power (Q_ESS) = P_total * tan(acos(PF))")
print(f"Q_ESS = {total_real_power:.2f} * tan(acos({power_factor}))")
print(f"Q_ESS = {total_reactive_power:.2f} MVAR")
print("\nThe required reactive power compensation from the ESS at Bus 9 is approximately 3.16 MVAR.")
