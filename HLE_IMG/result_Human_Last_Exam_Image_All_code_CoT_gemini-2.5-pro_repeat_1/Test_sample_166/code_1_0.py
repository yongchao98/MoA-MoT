import math

# Step 1: Define the given values for the loads and power factor.
P_L1 = 1.85  # MW
P_L2 = 1.70  # MW
P_L3 = 1.75  # MW
P_L4 = 1.90  # MW
P_L5 = 2.40  # MW
PF = 0.95    # Lagging power factor

print("This script calculates the required reactive power compensation from the ESS.")
print("-" * 60)

# Step 2: Calculate the total active power by summing all loads.
P_total = P_L1 + P_L2 + P_L3 + P_L4 + P_L5

print("1. Calculate the total active power (P_total):")
print(f"P_total = {P_L1} MW + {P_L2} MW + {P_L3} MW + {P_L4} MW + {P_L5} MW")
print(f"P_total = {P_total:.2f} MW")
print("-" * 60)

# Step 3: Calculate the required reactive power compensation.
# This is the total reactive power consumed by the loads, Q_total = P_total * tan(phi).
# First, find the angle phi from the power factor, phi = arccos(PF).
phi = math.acos(PF)
tan_phi = math.tan(phi)

# Then, calculate the reactive power.
Q_compensation = P_total * tan_phi

print("2. Calculate the required reactive power compensation (Q_compensation):")
print(f"The power factor (PF) is {PF}.")
print(f"The reactive power Q is calculated using the formula: Q = P_total * tan(arccos(PF))")
print("\nFinal Equation:")
print(f"Q_compensation = {P_total:.2f} MW * tan(arccos({PF})) = {P_total:.2f} MW * {tan_phi:.5f}")
print(f"Q_compensation = {Q_compensation:.3f} Mvar")
print("-" * 60)

print(f"The required reactive power compensation from the ESS is {Q_compensation:.3f} Mvar.")
<<<3.155>>>