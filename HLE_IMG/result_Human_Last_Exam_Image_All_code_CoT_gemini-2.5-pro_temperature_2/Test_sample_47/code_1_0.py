import numpy as np

# --- Define component values from the circuit diagram ---
R76 = 76
R8 = 8
R14 = 14
R11 = 11
R29 = 29
V_s = 41

# --- Set up the system of linear equations using KCL for nodal analysis ---
# The circuit has two unknown node voltages: VA (top wire) and VB (right junction).
# Equation for node A (sum of currents leaving = 0):
# (VA - V_s)/R76 + (VA - V_s)/R8 + (VA - VB)/R11 = 0
# Rearranged: VA*(1/R76 + 1/R8 + 1/R11) - VB*(1/R11) = V_s*(1/R76 + 1/R8)
#
# Equation for node B (sum of currents leaving = 0):
# (VB - VA)/R11 + (VB - V_s)/R14 + (VB - 0)/R29 = 0
# Rearranged: -VA*(1/R11) + VB*(1/R11 + 1/R14 + 1/R29) = V_s/R14

# --- Create the matrix (M) and vector (C) for the system M * v = C ---
# v = [VA, VB]
c11 = 1/R76 + 1/R8 + 1/R11
c12 = -1/R11
c21 = -1/R11
c22 = 1/R11 + 1/R14 + 1/R29

M = np.array([[c11, c12],
              [c21, c22]])

C = np.array([V_s * (1/R76 + 1/R8),
              V_s / R14])

# --- Solve for the node voltages VA and VB ---
try:
    voltages = np.linalg.solve(M, C)
    V_A = voltages[0]
    V_B = voltages[1]

    # --- Calculate the total current I leaving the source ---
    # I = current through R76 + current through R8 + current through R14
    I_R76 = (V_s - V_A) / R76
    I_R8 = (V_s - V_A) / R8
    I_R14 = (V_s - V_B) / R14
    I_total = I_R76 + I_R8 + I_R14

    # --- Print the results and the final equation with values ---
    print(f"Solving the circuit using nodal analysis gives the node voltages:")
    print(f"V_A (top wire voltage) = {V_A:.4f} V")
    print(f"V_B (junction voltage) = {V_B:.4f} V\n")

    print("The total current I is the sum of the currents leaving the 41V source.")
    print("The formula for I is: I = (V_s - V_A)/R76 + (V_s - V_A)/R8 + (V_s - V_B)/R14\n")

    print("Substituting the values into the equation:")
    print(f"I = ({V_s} - {V_A:.4f}) / {R76} + ({V_s} - {V_A:.4f}) / {R8} + ({V_s} - {V_B:.4f}) / {R14}")
    print(f"I = {I_R76:.4f} A + {I_R8:.4f} A + {I_R14:.4f} A")
    print(f"I = {I_total:.4f} A")

except np.linalg.LinAlgError:
    print("The system of equations cannot be solved. The matrix is singular.")
