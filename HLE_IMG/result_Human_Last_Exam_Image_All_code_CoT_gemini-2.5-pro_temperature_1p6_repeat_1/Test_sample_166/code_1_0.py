import math

# --- Given Data ---
# Load active powers in MW
P_L1 = 1.85
P_L2 = 1.7
P_L3 = 1.75
P_L4 = 1.9
P_L5 = 2.4

# Power factor for all loads (lagging)
PF = 0.95

# --- Step 1: Calculate the total active power load ---
P_total = P_L1 + P_L2 + P_L3 + P_L4 + P_L5

print("Step 1: Calculate the total active power load (P_total).")
print(f"P_total = {P_L1} + {P_L2} + {P_L3} + {P_L4} + {P_L5}")
print(f"P_total = {P_total:.2f} MW\n")

# --- Step 2: Calculate the required reactive power compensation ---
# The required reactive power (Q_ESS) is equal to the total reactive power
# consumed by the loads (Q_total).
# The formula is: Q_total = P_total * tan(acos(PF))
angle_phi = math.acos(PF)
tan_phi = math.tan(angle_phi)
Q_ESS = P_total * tan_phi

print("Step 2: Calculate the required reactive power compensation (Q_ESS).")
print("This is equal to the total reactive power demand of the loads.")
print(f"The calculation is based on the formula: Q_ESS = P_total * tan(acos(PF))")
print(f"Q_ESS = {P_total:.2f} MW * tan(acos({PF}))")
print(f"Q_ESS = {P_total:.2f} MW * {tan_phi:.4f}")
print(f"Q_ESS = {Q_ESS:.2f} MVAr\n")

print(f"The required reactive power compensation from the ESS is {Q_ESS:.2f} MVAr.")