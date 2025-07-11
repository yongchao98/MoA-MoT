import math

# Step 1: Extract Given Values
# Input power in Watts
P_in = 10e-3  # 10 mW
# Load resistance in Ohms
R_L = 2.7e3   # 2.7 KΩ
# Capacitor quality factor
Q_C = 150
# Operating frequency in MHz (for reading the plot)
f_MHz = 800

# Step 2: Determine Inductor Quality Factor (Q_L) from the plot
# From Figure (b), at a frequency of 800 MHz, we read the value on the y-axis
# for the dashed line (Quality factor). The value is approximately 98.
Q_L = 98

# Step 3: Calculate Power Transfer Efficiency (eta)
# Assuming "optimal power transfer" minimizes reflection losses, the main losses are
# from component quality factors. A plausible simplified model for the efficiency
# is the ratio of the component Q factors.
# eta = Q_inductor / Q_capacitor
eta = Q_L / Q_C

# Step 4: Calculate DC Output Power (P_DC)
# P_DC = P_in * eta
P_DC = P_in * eta

# Step 5: Calculate Voltage Across the Load (V_L)
# P_DC = V_L^2 / R_L  =>  V_L = sqrt(P_DC * R_L)
V_L = math.sqrt(P_DC * R_L)

# Output the results, showing the final equation and values
print("--- Calculation Steps ---")
print(f"1. Given values:")
print(f"   Input Power (P_in): {P_in:.3f} W")
print(f"   Load Resistance (R_L): {R_L:.0f} Ω")
print(f"   Capacitor Quality Factor (Q_C): {Q_C}")

print(f"\n2. Inductor Quality Factor (Q_L) from plot at {f_MHz} MHz: {Q_L}")

print(f"\n3. Power Transfer Efficiency (η) Calculation:")
print(f"   η = Q_L / Q_C = {Q_L} / {Q_C} = {eta:.4f}")

print(f"\n4. DC Output Power (P_DC) Calculation:")
print(f"   P_DC = P_in * η = {P_in:.3f} W * {eta:.4f} = {P_DC:.6f} W")

print(f"\n5. Final Voltage (V_L) Calculation:")
print(f"   V_L = sqrt(P_DC * R_L)")
print(f"   V_L = sqrt({P_DC:.6f} W * {R_L:.0f} Ω)")
print(f"   V_L = sqrt({P_DC * R_L:.4f}) V")
print(f"   V_L = {V_L:.2f} V")

print("\n--- Final Answer ---")
print(f"The calculated voltage across the load R_L is {V_L:.2f} V.")
print(f"Final equation: V_L = sqrt(({P_in}) * ({Q_L}/{Q_C}) * ({R_L}))")
print(f"Final calculation: V_L = sqrt({P_in} * {Q_L/Q_C:.4f} * {R_L}) = {V_L:.2f} V")
<<<4.20>>>