import math

# Step 1: Define the given parameters
P_in = 10e-3  # Input power in Watts (10 mW)
f = 0.8e9     # Operating frequency in Hertz (0.8 GHz)
Z0 = 50       # Characteristic impedance in Ohms
R_L = 2.7e3   # Load resistance in Ohms (2.7 kOhm)
Q_C = 150     # Quality factor of the capacitor (given)

# From the graph in Figure (b), at f = 800 MHz (0.8 GHz):
# The red dashed line (Quality factor) is approximately 95.
Q_L = 95      # Quality factor of the inductor from the plot

print("--- Parameters ---")
print(f"Input Power (P_in): {P_in * 1000} mW")
print(f"Operating Frequency (f): {f / 1e6} MHz")
print(f"Load Resistance (R_L): {R_L / 1000} kOhm")
print(f"Inductor Q-factor (Q_L) at {f / 1e6} MHz: {Q_L}")
print(f"Capacitor Q-factor (Q_C): {Q_C}")
print("-" * 20 + "\n")

# Step 2: Estimate the input resistance of the rectifier circuit
# For a half-wave rectifier, R_rect is approximately R_L / 2.
R_rect = R_L / 2
print("--- Calculations ---")
print(f"Step 2: Estimated rectifier input resistance (R_rect) = R_L / 2 = {R_rect} Ohms")

# Step 3: Calculate the loaded Q of the matching network
# The network matches R_rect to Z0.
# Q_load = sqrt(R_high / R_low - 1)
try:
    Q_load = math.sqrt(R_rect / Z0 - 1)
    print(f"Step 3: Calculated loaded Q (Q_load) = sqrt({R_rect:.1f} / {Z0} - 1) = {Q_load:.4f}")

    # Step 4: Calculate the efficiency of the matching network
    # eta_match is the product of efficiencies of its components.
    # We assume a two-component (L-C) matching network.
    # eta = (1 - Q_load/Q_L) * (1 - Q_load/Q_C)
    eta_match = (1 - Q_load / Q_L) * (1 - Q_load / Q_C)
    print(f"Step 4: Matching network efficiency (eta_match) = (1 - {Q_load:.4f}/{Q_L}) * (1 - {Q_load:.4f}/{Q_C}) = {eta_match:.4f}")

    # Step 5: Calculate the power delivered to the load
    # Assuming the rectifier is 100% efficient, P_load = P_in * eta_match
    P_load = P_in * eta_match
    print(f"Step 5: Power delivered to load (P_load) = {P_in} W * {eta_match:.4f} = {P_load * 1000:.4f} mW")

    # Step 6: Calculate the voltage across the load resistor
    V_L = math.sqrt(P_load * R_L)
    print(f"Step 6: Voltage across load (V_L) = sqrt(P_load * R_L) = sqrt({P_load:.6f} * {R_L}) = {V_L:.4f} V\n")
    
    # Final Answer Equation
    print("--- Final Equation ---")
    print(f"V_L = sqrt(P_in * (1 - Q_load/Q_L) * (1 - Q_load/Q_C) * R_L)")
    print(f"V_L = sqrt({P_in} * (1 - {Q_load:.3f}/{Q_L}) * (1 - {Q_load:.3f}/{Q_C}) * {R_L})")
    print(f"V_L = {V_L:.3f} V")

except ValueError:
    print("Error: Cannot calculate loaded Q. R_rect must be greater than Z0.")
