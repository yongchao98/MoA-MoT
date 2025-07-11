import math

# Step 1: Extract all given parameters
f = 0.8e9  # Operating frequency in Hz (0.8 GHz)
Pin = 10e-3  # Input power in Watts (10 mW)
QC = 150  # Quality factor of the capacitor
RL = 2.7e3  # Load resistance in Ohms (2.7 kOhm)
C1 = 0.4e-12 # Capacitance of C1 in Farads (0.4 pF)
L1 = 43e-9   # Inductance of L1 in Henries (43 nH)
# Step 2: Determine Inductor Quality Factor from the graph at 800 MHz
QL = 90  # Quality factor of the inductor from Figure (b) at 800 MHz

# Step 3 & 4: Calculate the total series loss resistance
# Angular frequency
omega = 2 * math.pi * f

# Reactance and series loss resistance for Inductor L1
XL1 = omega * L1
Rs_L1 = XL1 / QL

# Reactance and series loss resistance for Capacitor C1
XC1 = 1 / (omega * C1)
Rs_C1 = XC1 / QC # Or equivalently 1 / (omega * C1 * QC)

# Total series loss resistance
R_loss = Rs_L1 + Rs_C1

# Step 5: Calculate the power transfer efficiency based on the simplified loss model
# This model assumes power is divided between the series loss resistance and the load resistance
eta = RL / (RL + R_loss)

# Step 6: Calculate the power delivered to the load
P_load = Pin * eta

# Step 7: Calculate the voltage across the load
V_L = math.sqrt(P_load * RL)

# Print the results step-by-step
print(f"--- Calculation Steps ---")
print(f"Operating Frequency (f): {f/1e6} MHz")
print(f"Input Power (Pin): {Pin*1000} mW")
print(f"Inductor L1: {L1*1e9} nH, Quality Factor (QL): {QL}")
print(f"Capacitor C1: {C1*1e12} pF, Quality Factor (QC): {QC}")
print(f"Load Resistance (RL): {RL/1000} kOhms")
print(f"\n1. Calculate series loss resistance:")
print(f"   - Loss resistance of L1 (Rs_L1): {Rs_L1:.4f} Ohms")
print(f"   - Loss resistance of C1 (Rs_C1): {Rs_C1:.4f} Ohms")
print(f"   - Total Loss Resistance (R_loss): {R_loss:.4f} Ohms")

print(f"\n2. Calculate efficiency and power at the load:")
print(f"   - Power Transfer Efficiency (eta): {eta:.4f}")
print(f"   - Power at Load (P_load = Pin * eta): {P_load*1000:.4f} mW")

print(f"\n3. Calculate the final voltage across the load RL:")
print(f"P_load = {P_load:.6f} W")
print(f"R_L = {RL} Ohms")
print(f"V_L = sqrt(P_load * R_L) = sqrt({P_load:.6f} * {RL}) = {V_L:.4f} V")
print(f"\nFinal calculated voltage across the load is {V_L:.4f} V.")

# Final answer in the required format
# <<<5.1907>>>