import numpy as np

# --- Define component values ---
R76 = 76.0
R8 = 8.0
R14 = 14.0
R11 = 11.0
R29 = 29.0
V_source = 41.0

# --- Nodal Analysis ---
# The circuit diagram interpretation leads to a system of equations for node potentials.
# Nodes are defined as:
# V_A = 41V (top rail)
# V_B = 0V (bottom rail/ground)
# V_C = node between R76 and R8
# V_D = node between R8, R14, and the wire with current I
# V_E = node between R14, R11, R29
# The wire with current I connects ground (V_B) to node D.
# Assuming an ideal wire, V_D = V_B = 0V.

# 1. Solve for V_C
# KCL at Node C: (V_C - V_A)/R76 + (V_C - V_D)/R8 = 0
# With V_A = 41 and V_D = 0: (V_C - 41)/76 + V_C/8 = 0
# V_C * (1/76 + 1/8) = 41/76
V_C = (V_source / R76) / (1/R76 + 1/R8)

# 2. Solve for V_E
# KCL at Node E: (V_E - V_A)/R11 + (V_E - V_D)/R14 + (V_E - V_B)/R29 = 0
# With V_A = 41, V_D = 0, V_B = 0: (V_E - 41)/11 + V_E/14 + V_E/29 = 0
# V_E * (1/11 + 1/14 + 1/29) = 41/11
V_E = (V_source / R11) / (1/R11 + 1/R14 + 1/R29)

# 3. Solve for I
# I is the current entering Node D from ground.
# KCL at Node D: I_from_C + I_from_E + I = 0
# (V_C - V_D)/R8 + (V_E - V_D)/R14 + I = 0
# With V_D = 0: V_C/R8 + V_E/R14 + I = 0
# I = -V_C/R8 - V_E/R14
I = -V_C / R8 - V_E / R14

# --- Print the results step-by-step ---
print("Step 1: Calculate the potential at Node C (V_C)")
print(f"V_C = ({V_source} / {R76}) / (1/{R76} + 1/{R8}) = {V_C:.3f} V")
print("\nStep 2: Calculate the potential at Node E (V_E)")
print(f"V_E = ({V_source} / {R11}) / (1/{R11} + 1/{R14} + 1/{R29}) = {V_E:.3f} V")
print("\nStep 3: Calculate the current I using KCL at Node D")
print("I = -V_C / R8 - V_E / R14")
print(f"I = -({V_C:.3f} V) / {R8} \u03A9 - ({V_E:.3f} V) / {R14} \u03A9")
print(f"I = {-V_C/R8:.3f} A - {V_E/R14:.3f} A")
print(f"I = {I:.3f} A")

# Final answer in the required format
# print(f"\nFinal Answer: {I:.1f}")