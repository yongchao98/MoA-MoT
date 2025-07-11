import numpy as np

# Step 1: Define the component values from the circuit diagram.
R_76 = 76  # ohms
R_8 = 8    # ohms
R_14 = 14  # ohms
R_11 = 11  # ohms
R_29 = 29  # ohms
V_source = 41 # Volts

# Step 2: Set up the system of linear equations (Ax = b) based on Nodal Analysis.
# Let's define the unknown node voltages:
# V1: Potential of the top-most wire.
# V2: Potential of the wire where current 'I' flows.
# VR: Potential of the node connecting R_14, R_11, and R_29.
#
# Known voltages:
# V_source_line = 41V (wire connected to the positive terminal of the battery)
# V_ground = 0V (bottom wire)
#
# KCL Equations:
# Node V1: (V1 - V2)/R_76 + (V1 - VR)/R_11 = 0
# Node V2: (V2 - V1)/R_76 + (V2 - V_source_line)/R_8 + (V2 - VR)/R_14 = 0
# Node VR: (VR - V1)/R_11 + (VR - V2)/R_14 + (VR - V_ground)/R_29 = 0
#
# Rearranging into Ax = b form:
# V1*(1/76+1/11) - V2*(1/76) - VR*(1/11) = 0
# -V1*(1/76) + V2*(1/76+1/8+1/14) - VR*(1/14) = 41/8
# -V1*(1/11) - V2*(1/14) + VR*(1/11+1/14+1/29) = 0

# Matrix A represents the coefficients of the unknown voltages.
A = np.array([
    [1/R_76 + 1/R_11, -1/R_76, -1/R_11],
    [-1/R_76, 1/R_76 + 1/R_8 + 1/R_14, -1/R_14],
    [-1/R_11, -1/R_14, 1/R_11 + 1/R_14 + 1/R_29]
])

# Vector b represents the constants on the right side of the equations.
b = np.array([0, V_source/R_8, 0])

# Step 3: Solve for the unknown node voltages x = [V1, V2, VR].
try:
    voltages = np.linalg.solve(A, b)
    V1, V2, VR = voltages[0], voltages[1], voltages[2]

    # Step 4: Calculate the current I flowing through the 14-ohm resistor.
    # The arrow for I indicates current from node V2 to node VR.
    # I = (V2 - VR) / R_14
    current_I = (V2 - VR) / R_14

    # Step 5: Print the results, showing the final equation with values.
    print(f"Calculated Node Voltages:")
    print(f"V_left (V2) = {V2:.4f} V")
    print(f"V_right (VR) = {VR:.4f} V")
    print("\nCalculating the current I:")
    print(f"I = (V_left - V_right) / R_14")
    # Here is the final equation with all numbers plugged in
    print(f"I = ({V2:.4f} - {VR:.4f}) / {R_14}")
    print(f"The value of the current I is: {current_I:.4f} A")

except np.linalg.LinAlgError:
    print("Could not solve the system of equations. The matrix is singular.")
