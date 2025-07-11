import math

# Given values
P_in = 10e-3  # Input power in Watts (10 mW)
R_L = 2.7e3   # Load resistance in Ohms (2.7 kOhm)

# The problem states the circuit is designed for "optimal power transfer" and to "maximize power transfer efficiency".
# For a network with internal losses (due to component Q-factors), the condition for maximum power transfer
# occurs when the power delivered to the load equals the power dissipated in the network's lossy components.
# This results in a theoretical efficiency of 50%.
eta_total = 0.5

# Calculate the DC power delivered to the load
P_L = P_in * eta_total

# Calculate the voltage across the load
V_L = math.sqrt(P_L * R_L)

# Print the final equation with all the numbers
print(f"The calculation for the load voltage V_L is based on the formula V_L = sqrt(P_in * eta_total * R_L)")
print(f"V_L = sqrt({P_in} W * {eta_total} * {R_L} Ohms)")
print(f"V_L = sqrt({P_L} W * {R_L} Ohms)")
print(f"V_L = sqrt({P_L * R_L}) V")
print(f"The calculated voltage across the load is: {V_L:.4f} V")

# Final answer in the required format
print(f"<<<{V_L:.4f}>>>")