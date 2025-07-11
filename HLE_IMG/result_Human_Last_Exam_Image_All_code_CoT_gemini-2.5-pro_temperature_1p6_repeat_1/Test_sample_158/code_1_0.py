import math

# Step 1 & 2: Define constants from the problem description and the graph.
# Input power in Watts
P_in = 10e-3  # 10 mW
# Load resistance in Ohms
R_L = 2.7e3   # 2.7 kOhm
# From Figure (b), at the operating frequency of 800 MHz (0.8 GHz),
# the quality factor of the inductor (dashed red line) is approximately 90.
Q_L = 90
# The quality factor of the capacitor is given.
Q_C = 150

# Step 3 & 4: Model the losses and calculate the total efficiency.
# The problem states the circuit is optimized for power transfer, implying reflection losses are minimized (assumed to be 0).
# The primary losses are from the finite quality factors of the components.
# We will calculate a total efficiency based on these Q-factors, assuming ideal diode conversion.
# Loss tangent (dissipation factor) d = 1/Q.
# A simplified model for total loss is the sum of individual loss factors.
d_L = 1 / Q_L
d_C = 1 / Q_C
total_loss_factor = d_L + d_C
eta_total = 1 - total_loss_factor

# Step 5: Calculate the output power delivered to the load.
P_L = P_in * eta_total

# Step 6: Calculate the voltage across the load resistor.
V_L = math.sqrt(P_L * R_L)

# Step 7: Print the final output as requested.
# The final equation is V_L = sqrt(P_in * eta_total * R_L)
print(f"Based on the given data and a simplified loss model, we calculate the voltage across the load R_L.")
print(f"Inductor Quality Factor (Q_L) at 800 MHz: {Q_L}")
print(f"Capacitor Quality Factor (Q_C): {Q_C}")
print(f"Total Efficiency (eta) = 1 - (1/Q_L + 1/Q_C) = {eta_total:.4f}")
print(f"Output Power (P_L) = P_in * eta = {P_L*1000:.4f} mW")
print("\nFinal Voltage Calculation:")
# Print the final equation with all the numbers
print(f"V_L = sqrt(P_in * eta * R_L)")
print(f"V_L = sqrt({P_in:.3f} W * {eta_total:.4f} * {R_L:.1f} Ω)")
print(f"V_L = {V_L:.4f} V")

# The final answer in the required format.
# Round the final voltage to two decimal places for the answer.
final_answer = round(V_L, 2)
# print(f"\nFinal Answer: {final_answer}")
print(f'\n<<<sqrt({P_in} * {eta_total} * {R_L}) = {final_answer}>>>')