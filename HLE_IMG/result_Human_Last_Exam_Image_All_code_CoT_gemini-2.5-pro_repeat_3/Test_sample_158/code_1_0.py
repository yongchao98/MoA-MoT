import math

# Step 1 & 2: Define the given parameters
# Input power in Watts
P_in = 10e-3  # 10 mW
# Load resistance in Ohms
R_L = 2.7e3  # 2.7 kΩ
# Capacitor quality factor
Q_C = 150
# Inductor quality factor from the graph at 800 MHz
Q_L = 95

# Step 3: Model and calculate the overall power conversion efficiency (PCE)
# Based on the assumption that efficiency is determined by the ratio of component Q-factors.
# η_PCE = 1 / (1 + (Power_Loss_Ratio))
# Assuming Power_Loss_Ratio is proportional to Q_L / Q_C
eta_PCE = 1 / (1 + Q_L / Q_C)

# Step 4: Calculate the power delivered to the load
P_L = P_in * eta_PCE

# Step 5: Calculate the voltage across the load
V_L = math.sqrt(P_L * R_L)

# Output the final equation with numerical values and the result
print("Based on the simplified efficiency model, the calculation for the load voltage V_L is:")
print(f"V_L = sqrt((1 / (1 + {Q_L} / {Q_C})) * {P_in:.3f} W * {R_L:.0f} Ω)")
print(f"V_L = sqrt({eta_PCE:.4f} * {P_in:.3f} W * {R_L:.0f} Ω)")
print(f"V_L = sqrt({P_L:.6f} W * {R_L:.0f} Ω)")
print(f"V_L = {V_L:.2f} V")

# Final answer in the required format
# The final result is the voltage V_L
final_answer = V_L
print(f'<<<{final_answer:.2f}>>>')