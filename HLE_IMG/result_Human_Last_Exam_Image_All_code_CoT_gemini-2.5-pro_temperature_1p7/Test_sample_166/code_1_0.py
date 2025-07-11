import math

# Step 1: Define the real power of each load in MW from the diagram.
P_L1 = 1.85
P_L2 = 1.7
P_L3 = 1.75
P_L4 = 1.9
P_L5 = 2.4

# Step 2: Calculate the total real power demand (P_total) by summing the individual loads.
P_total = P_L1 + P_L2 + P_L3 + P_L4 + P_L5

# Step 3: Define the given power factor (PF).
PF = 0.95

# Step 4: Calculate the total reactive power demand (Q_total).
# The relationship is Q_total = P_total * tan(arccos(PF)).
# First, calculate the value of tan(arccos(PF)).
phi = math.acos(PF)
tan_phi = math.tan(phi)
# Then, calculate the total reactive power.
Q_total = P_total * tan_phi

# The required reactive power compensation from the ESS is equal to the total reactive power demand.
Q_compensation = Q_total

print("Step 1: Calculate the total real power (P_total) demand.")
print(f"P_total = {P_L1} + {P_L2} + {P_L3} + {P_L4} + {P_L5} = {P_total:.2f} MW")
print("-" * 30)

print("Step 2: Use the power factor (PF) to find the total reactive power (Q_total).")
print(f"The formula is: Q_total = P_total * tan(arccos(PF))")
print(f"Q_total = {P_total:.2f} * tan(arccos({PF}))")
print(f"Q_total = {P_total:.2f} * {tan_phi:.4f} = {Q_total:.2f} MVAR")
print("-" * 30)

print("Step 3: Determine the required compensation from the ESS.")
print("To correct the voltage drop, the ESS must compensate for the total reactive power demand of the system.")
print(f"Required Reactive Power Compensation = {Q_compensation:.2f} MVAR")
print("-" * 30)

# Final answer in the equation format
print("Final Calculation:")
print(f"Q_compensation = ({P_L1} + {P_L2} + {P_L3} + {P_L4} + {P_L5}) * tan(arccos({PF})) = {Q_compensation:.2f} MVAR")