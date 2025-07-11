import numpy as np

# --- 1. Define component values ---
R1 = 76.0  # Ohms
R2 = 8.0   # Ohms
R3 = 14.0  # Ohms
R4 = 11.0  # Ohms
R5 = 29.0  # Ohms
V_source = 41.0 # Volts

# --- 2. Set up the system of linear equations for Nodal Analysis (Ax = b) ---
# For nodes V_L and V_R

# KCL at V_L: V_L*(1/R1 + 1/R2 + 1/R3) - V_R*(1/R3) = V_source*(1/R1 + 1/R2)
# KCL at V_R: -V_L*(1/R3) + V_R*(1/R4 + 1/R5 + 1/R3) = V_source*(1/R4)

# Matrix A
A = np.array([
    [(1/R1 + 1/R2 + 1/R3), -1/R3],
    [-1/R3, (1/R4 + 1/R5 + 1/R3)]
])

# Vector b
b = np.array([
    V_source * (1/R1 + 1/R2),
    V_source / R4
])

# --- 3. Solve for the node voltages V_L and V_R ---
try:
    voltages = np.linalg.solve(A, b)
    V_L = voltages[0]
    V_R = voltages[1]

    # --- 4. Calculate the current I ---
    I = (V_L - V_R) / R3

    # --- 5. Print the results ---
    print(f"Based on the analysis, the node voltages are:")
    print(f"V_L = {V_L:.3f} V")
    print(f"V_R = {V_R:.3f} V")
    print("\nThe current I is the current through the 14 Ohm resistor, from left to right.")
    print(f"I = (V_L - V_R) / R14")
    print(f"I = ({V_L:.3f} V - {V_R:.3f} V) / {R3} Ohm")
    print(f"I = {I:.3f} A")

except np.linalg.LinAlgError:
    print("The system of equations could not be solved. The circuit might be ill-defined.")
