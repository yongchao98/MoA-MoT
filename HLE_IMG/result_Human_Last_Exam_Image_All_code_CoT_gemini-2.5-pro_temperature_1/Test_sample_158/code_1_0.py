import math

# Step 1: Define the given values based on the problem description and images.
P_in = 10e-3  # Input power in Watts (10 mW)
R_L = 2.7e3   # Load resistance in Ohms (2.7 kOhm)
Q_C = 150     # Quality factor of the capacitor (given)

# Step 2: Determine the inductor quality factor from the graph at 0.8 GHz (800 MHz).
# From the graph, the dashed red line (Quality factor) is at 90 at 800 MHz.
Q_L = 90

# Step 3: Calculate the overall power transfer efficiency (eta) based on component losses.
# The model assumes efficiency is limited by the Q-factors of the inductor and capacitor.
eta = (1 - 1/Q_L) * (1 - 1/Q_C)

# Step 4: Calculate the output power delivered to the load.
P_out = P_in * eta

# Step 5: Calculate the voltage across the load resistor.
V_L = math.sqrt(P_out * R_L)

# Print the final equation with all numerical values substituted, as requested.
print("Calculating the voltage across the load (V_L):")
print(f"V_L = sqrt(P_in * (1 - 1/Q_L) * (1 - 1/Q_C) * R_L)")
print(f"V_L = sqrt({P_in} * (1 - 1/{Q_L}) * (1 - 1/{Q_C}) * {R_L})")

# Print the final calculated result.
print(f"\nThe calculated voltage across the load is: {V_L:.2f} V")

# The final answer in the required format
# print(f'<<<{V_L:.2f}>>>')