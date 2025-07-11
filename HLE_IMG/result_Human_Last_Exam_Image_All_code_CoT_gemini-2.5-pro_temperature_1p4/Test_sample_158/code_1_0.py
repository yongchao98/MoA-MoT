import math

# --- Given Parameters ---
# Input RF Power in Watts
P_in = 10e-3  # 10 mW
# Load Resistance in Ohms
R_L = 2.7e3   # 2.7 kOhm

# --- Assumption ---
# The problem states the circuit is designed for "optimal power transfer".
# For a resonant system, this implies critical coupling, where the theoretical
# power transfer efficiency is 50%.
eta = 0.5

# --- Calculation ---
# 1. Calculate the DC power delivered to the load
P_L = P_in * eta

# 2. Calculate the voltage across the load
V_L = math.sqrt(P_L * R_L)

# --- Output the results step-by-step ---
print("Step 1: Define input values and assumption")
print(f"Input Power (P_in) = {P_in} W")
print(f"Load Resistance (R_L) = {R_L} Ohms")
print(f"Assumed Efficiency for Optimal Power Transfer (eta) = {eta}")
print("-" * 30)

print("Step 2: Calculate Power on Load (P_L)")
print(f"P_L = P_in * eta")
print(f"P_L = {P_in} W * {eta} = {P_L} W")
print("-" * 30)

print("Step 3: Calculate Voltage on Load (V_L)")
print(f"V_L = sqrt(P_L * R_L)")
print(f"V_L = sqrt({P_L} * {R_L})")
print(f"V_L = {V_L:.4f} V")

# Final Answer in the requested format
# print(f"<<<{V_L:.4f}>>>")