import math

# Step 1: Define the given parameters
P_in = 10e-3  # Input power in Watts (10 mW)
f = 0.8e9     # Operating frequency in Hz (0.8 GHz)
R_L = 2.7e3   # Load resistance in Ohms (2.7 KΩ)
Q_C = 150     # Quality factor of the capacitor

# Step 2: Determine the inductor's quality factor from the graph
# At f = 800 MHz, the red dashed line (Quality factor) in Figure (b) is at 90.
Q_L = 90

print(f"Given parameters:")
print(f"Input Power (Pin): {P_in * 1000} mW")
print(f"Load Resistance (RL): {R_L / 1000} KΩ")
print(f"Capacitor Quality Factor (QC): {Q_C}")
print(f"Inductor Quality Factor (QL) at {f/1e6} MHz from graph: {Q_L}\n")

# Step 3: Calculate the total quality factor (Q_total)
# The loss effects of the inductor and capacitor combine.
# 1/Q_total = 1/Q_L + 1/Q_C
Q_total_inv = (1 / Q_L) + (1 / Q_C)
Q_total = 1 / Q_total_inv
print(f"Calculating the total quality factor (Q_total):")
print(f"1 / Q_total = 1 / {Q_L} + 1 / {Q_C}")
print(f"Q_total = {Q_total:.2f}\n")

# Step 4: Assume a model for Power Conversion Efficiency (PCE)
# A plausible simplification is that the PCE is related to the total Q factor.
# We'll assume PCE = Q_total / 100.
eta_PCE = Q_total / 100
print(f"Calculating the Power Conversion Efficiency (PCE or η):")
print(f"Assuming η = Q_total / 100")
print(f"η = {Q_total:.2f} / 100 = {eta_PCE:.4f} or {eta_PCE*100:.2f}%\n")

# Step 5: Calculate the output power (P_out)
# P_out = P_in * η_PCE
P_out = P_in * eta_PCE
print(f"Calculating the output power (P_out):")
print(f"P_out = Pin * η = {P_in:.3f} W * {eta_PCE:.4f}")
print(f"P_out = {P_out:.6f} W\n")

# Step 6: Calculate the voltage across the load (V_L)
# P_out = V_L^2 / R_L  =>  V_L = sqrt(P_out * R_L)
V_L = math.sqrt(P_out * R_L)
print(f"Calculating the voltage across the load (V_L):")
print(f"V_L = sqrt(P_out * R_L)")
print(f"V_L = sqrt({P_out:.6f} W * {R_L:.0f} Ω)")
print(f"V_L = {V_L:.4f} V")

# Final answer
final_voltage = V_L
print(f"<<<{final_voltage:.2f}>>>")