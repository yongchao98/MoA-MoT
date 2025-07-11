import numpy as np

# --- 1. Define resistor values and source voltage ---
R76 = 76
R8 = 8
R14 = 14
R11 = 11
R29 = 29
V_s = 41

# --- 2. Set up the system of linear equations Ax = B for V1 and V2 ---
# Equation 1 for node V1: V1*(1/R8 + 1/R76 + 1/R14) - V2*(1/R14) = V_s/R8
# Equation 2 for node V2: -V1*(1/R14) + V2*(1/R29 + 1/R11 + 1/R14) = V_s/R29

A = np.array([
    [1/R8 + 1/R76 + 1/R14, -1/R14],
    [-1/R14, 1/R29 + 1/R11 + 1/R14]
])

B = np.array([
    V_s/R8,
    V_s/R29
])

# --- 3. Solve for V = [V1, V2] ---
try:
    V = np.linalg.solve(A, B)
    V1 = V[0]
    V2 = V[1]
    
    # --- 4. Calculate the current I ---
    # I is the current flowing into the second part of the circuit, which is the current through R29
    I = (V_s - V2) / R29

    # --- 5. Print the results ---
    print(f"Solving the KCL equations for node voltages V1 and V2.")
    print(f"V1 = {V1:.4f} V")
    print(f"V2 = {V2:.4f} V")
    print("\nThe current I is the current flowing from the first branch takeoff to the second.")
    print("This corresponds to the current passing through the 29 Ohm resistor.")
    print(f"\nI = (V_source - V2) / R29")
    print(f"I = ({V_s} - {V2:.4f}) / {R29}")
    print(f"I = {V_s - V2:.4f} / {R29}")
    print(f"I = {I:.4f} A")

except np.linalg.LinAlgError:
    print("The system of equations is singular and cannot be solved.")
