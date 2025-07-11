import numpy as np

# --- Given values ---
V_source = 41.0  # Volts
R1 = 76.0  # Ohms (Resistor between 41V and Node A)
R2 = 8.0   # Ohms (Resistor for current I, between 41V and Node A)
R3 = 14.0  # Ohms (Resistor between Node A and Node B)
R4 = 11.0  # Ohms (Resistor between 41V and Node B)
R5 = 29.0  # Ohms (Resistor between Node B and Ground)

# --- Nodal Analysis Equations ---
# Equation for Node A: Va*(1/R1 + 1/R2 + 1/R3) - Vb*(1/R3) = V_source/R1 + V_source/R2
# Equation for Node B: -Va*(1/R3) + Vb*(1/R3 + 1/R4 + 1/R5) = V_source/R4

# --- Setup the system of linear equations Ax = B ---
# A is the matrix of coefficients for Va and Vb
A = np.array([
    [1/R1 + 1/R2 + 1/R3, -1/R3],
    [-1/R3, 1/R3 + 1/R4 + 1/R5]
])

# B is the vector of constants
B = np.array([
    V_source/R1 + V_source/R2,
    V_source/R4
])

# --- Solve for potentials [Va, Vb] ---
potentials = np.linalg.solve(A, B)
Va = potentials[0]
Vb = potentials[1]

# --- Calculate the current I ---
# I is the current through R2
I = (V_source - Va) / R2

# --- Print the results ---
print(f"The potential at Node A (Va) is: {Va:.4f} V")
print(f"The potential at Node B (Vb) is: {Vb:.4f} V")
print("\nThe current I is calculated using Ohm's Law for the 8 Ohm resistor:")
print(f"I = (V_source - Va) / R2")
print(f"I = ({V_source} V - {Va:.4f} V) / {R2} Ohm")
print(f"I = {I:.4f} A")