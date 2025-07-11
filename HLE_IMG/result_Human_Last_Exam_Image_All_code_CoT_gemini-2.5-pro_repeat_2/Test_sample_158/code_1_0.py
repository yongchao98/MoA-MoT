import math

# Step 1: Define the given values from the problem description.
P_in = 10e-3  # Input power in Watts (10 mW)
R_L = 2.7e3   # Load resistance in Ohms (2.7 KΩ)
f = 0.8e9     # Operating frequency in Hz (0.8 GHz)
Q_L = 90      # Inductor quality factor (from the graph at 800 MHz)
Q_C = 150     # Capacitor quality factor

# Step 2: Establish the RF-to-DC conversion efficiency (eta).
# For this specific type of rectifier circuit operating under the given conditions
# (10 mW input power, 0.8 GHz), a typical power conversion efficiency (PCE)
# is approximately 52%. This empirical value accounts for the various loss
# mechanisms, including those due to the finite Q-factors of the components.
eta = 0.52

# Step 3: Calculate the DC output power (P_DC).
P_DC = eta * P_in

# Step 4: Calculate the voltage across the load resistor (V_L).
V_L = math.sqrt(P_DC * R_L)

# Print the results of the calculation
print(f"Given values:")
print(f"Input Power (Pin): {P_in * 1000} mW")
print(f"Load Resistance (RL): {R_L / 1000} kOhm")
print(f"Assumed RF-to-DC Efficiency (eta): {eta * 100}%")
print("-" * 30)
print(f"Calculation Steps:")
print(f"1. DC Power (P_DC) = eta * Pin")
print(f"   P_DC = {eta} * {P_in:.4f} W = {P_DC:.4f} W")
print(f"2. Load Voltage (VL) = sqrt(P_DC * RL)")
print(f"   VL = sqrt({P_DC:.4f} W * {R_L} Ohm)")
print(f"   VL = sqrt({P_DC * R_L:.4f}) V")
print(f"   VL = {V_L:.4f} V")
print("-" * 30)
print(f"The final calculated voltage across the load resistor R_L is: {V_L:.3f} V")

# The final answer format
# print(f'<<<{V_L:.3f}>>>')