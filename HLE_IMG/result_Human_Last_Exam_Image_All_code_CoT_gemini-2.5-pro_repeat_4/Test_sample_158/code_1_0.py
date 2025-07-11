import math

# Step 1 & 2: Define constants and parameters from the problem and graph
P_in = 10e-3      # Input power in Watts (10 mW)
f = 0.8e9         # Operating frequency in Hertz (0.8 GHz)
R_L = 2.7e3       # Load resistance in Ohms (2.7 kOhm)
Z_0 = 50          # Characteristic impedance in Ohms
Q_C = 150         # Quality factor of capacitors
Q_L = 95          # Quality factor of inductors from graph at 800 MHz

# Component values from the table
L1 = 43e-9        # Inductance of L1 in Henrys (43 nH)
C1 = 0.4e-12      # Capacitance of C1 in Farads (0.4 pF)
L2 = 39e-9        # Inductance of L2 in Henrys (39 nH)
C2 = 0.2e-12      # Capacitance of C2 in Farads (0.2 pF)

# Step 3: Calculate equivalent series resistance (Rs) for each component
omega = 2 * math.pi * f

# Loss resistance for inductors
Rs_L1 = (omega * L1) / Q_L
Rs_L2 = (omega * L2) / Q_L

# Loss resistance for capacitors
Rs_C1 = 1 / (omega * C1 * Q_C)
Rs_C2 = 1 / (omega * C2 * Q_C)

# Step 4: Calculate the equivalent loss resistance of the network
Rs_branch1 = Rs_L1 + Rs_C1
Rs_branch2 = Rs_L2 + Rs_C2
Rs_eq = (Rs_branch1 + Rs_branch2) / 2

# Step 5: Determine the equivalent load resistance from the matching condition
# Z_in_rect = R_load_eq + Rs_eq = Z_0
R_load_eq = Z_0 - Rs_eq

# Step 6: Calculate the power transfer efficiency (eta)
eta = R_load_eq / Z_0

# Step 7: Calculate the output power (P_out)
P_out = P_in * eta

# Step 8: Calculate the voltage across the load (V_L)
V_L = math.sqrt(P_out * R_L)

# Step 9: Print the final calculation steps and the result
print("Final Calculation Steps:")
print(f"1. Calculated Power Transfer Efficiency (eta) = {eta:.4f}")
print(f"2. Given Input Power (Pin) = {P_in:.3f} W")
print(f"3. Given Load Resistance (RL) = {R_L:.0f} Ohms")
print("\nFinal Equation:")
print(f"V_L = sqrt(eta * Pin * RL)")
print(f"V_L = sqrt({eta:.4f} * {P_in:.3f} * {R_L:.0f})")
print(f"The calculated voltage across the load is {V_L:.4f} V.")

<<<4.8149>>>