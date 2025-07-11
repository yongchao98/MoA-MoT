import numpy as np

# --- Define the known values from the circuit diagram ---
V_source = 41  # Voltage of the source in Volts
R76 = 76       # Resistance in Ohms
R8 = 8         # Resistance in Ohms
R14 = 14       # Resistance in Ohms
R11 = 11       # Resistance in Ohms
R29 = 29       # Resistance in Ohms

# --- Set up the system of linear equations (A * x = b) for nodal analysis ---
# where x = [V_B, V_D] are the unknown node voltages.

# The coefficient matrix A is derived from the KCL equations:
# Eq1: V_B*(1/R76+1/R8+1/R14) - V_D*(1/R14) = V_source/R76 + V_source/R8
# Eq2: -V_B*(1/R14) + V_D*(1/R11+1/R29+1/R14) = V_source/R11
A = np.array([
    [(1/R76 + 1/R8 + 1/R14), -1/R14],
    [-1/R14, (1/R11 + 1/R29 + 1/R14)]
])

# The constant vector b is also derived from the KCL equations:
b = np.array([
    V_source/R76 + V_source/R8,
    V_source/R11
])

# --- Solve for the unknown node voltages V_B and V_D ---
try:
    node_voltages = np.linalg.solve(A, b)
    V_B = node_voltages[0]
    V_D = node_voltages[1]

    # --- Calculate the current I using Ohm's Law ---
    # I flows from node B to node D through R14.
    I = (V_B - V_D) / R14

    # --- Print the results ---
    print(f"Solved Node Voltages:")
    print(f"Voltage at node B (V_B): {V_B:.2f} V")
    print(f"Voltage at node D (V_D): {V_D:.2f} V")
    print("\nCalculating the current I:")
    # Using f-string to format the final equation with calculated values
    print(f"I = (V_B - V_D) / R14 = ({V_B:.2f} V - {V_D:.2f} V) / {R14} Ω = {I:.3f} A")

except np.linalg.LinAlgError:
    print("The system of equations cannot be solved. The matrix is singular.")
