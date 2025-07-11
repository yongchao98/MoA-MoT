import numpy as np

# Step 1: Define circuit parameters from the schematic
V_s = 41.0  # Source Voltage (V)
R_76 = 76.0 # Resistor (Ohms)
R_8 = 8.0   # Resistor (Ohms)
R_14 = 14.0 # Resistor (Ohms)
R_11 = 11.0 # Resistor (Ohms)
R_29 = 29.0 # Resistor (Ohms)
V_gnd = 0.0 # Ground Voltage (V)

print("Circuit Analysis Steps:")
print("1. Define circuit topology and node potentials.")
print("   - Source Vs = 41V, Ground = 0V.")
print("   - R_76 is between Vs and Node B.")
print("   - R_8 is between Vs and Node A.")
print("   - R_14 is between Node A and Node B.")
print("   - R_11 is between Node B and Node C.")
print("   - R_29 is between Node C and Ground.")
print("   - The wire with current 'I' connects Vs to Node C, so Vc = 41V.")
print("   - The unknown potentials are Va and Vb.\n")

print("2. Set up and solve KCL equations for Va and Vb.")
# We can represent the system of equations as M*v = b, where v = [Va, Vb]'
# KCL at Node A: (Va-Vs)/R_8 + (Va-Vb)/R_14 = 0
# KCL at Node B: (Vb-Vs)/R_76 + (Vb-Va)/R_14 + (Vb-Vs)/R_11 = 0 (since Vc=Vs)

# Conductance matrix M
G_76 = 1/R_76
G_8 = 1/R_8
G_14 = 1/R_14
G_11 = 1/R_11

M = np.array([
    [G_8 + G_14, -G_14],
    [-G_14, G_76 + G_14 + G_11]
])

# Current vector b
b = np.array([
    V_s * G_8,
    V_s * G_76 + V_s * G_11
])

# Solve for v = [Va, Vb]
try:
    v = np.linalg.solve(M, b)
    Va = v[0]
    Vb = v[1]
except np.linalg.LinAlgError:
    print("Error: The matrix is singular, cannot solve.")
    Va, Vb = 0, 0
    
Vc = V_s

print(f"   Solved potentials: Va = {Va:.2f}V, Vb = {Vb:.2f}V, Vc = {Vc:.2f}V\n")


print("3. Calculate the current I using KCL at Node C.")
print("   The current 'I' flowing into Node C equals the sum of currents leaving it.")
# Current from C to B through R_11
I_R11 = (Vc - Vb) / R_11
# Current from C to Ground through R_29
I_R29 = (Vc - V_gnd) / R_29
# Total current I
I = I_R11 + I_R29

print(f"   I = (Vc - Vb)/R_11 + (Vc - Vgnd)/R_29")
print(f"   I = ({Vc:.2f} V - {Vb:.2f} V) / {R_11} Ω + ({Vc:.2f} V - {V_gnd:.2f} V) / {R_29} Ω")
print(f"   I = {I_R11:.4f} A + {I_R29:.4f} A")
print(f"   I = {I:.4f} A\n")

print("Final Answer:")
print(f"The value of the current I is {I:.2f} A.")
