import numpy as np

# Define the given values from the circuit diagram
R76 = 76.0
R8 = 8.0
R14 = 14.0
R11 = 11.0
R29 = 29.0
V_source = 41.0

print("Step 1: Simplify the parallel resistors on the left side.")
# The 76 Ohm and 8 Ohm resistors are in parallel. Calculate their equivalent resistance, Rp.
Rp = (R76 * R8) / (R76 + R8)
print(f"The equivalent resistance of the {R76} Ohm and {R8} Ohm resistors in parallel is:")
print(f"Rp = ({R76} * {R8}) / ({R76} + {R8}) = {Rp:.4f} Ohms")
print("-" * 30)

print("Step 2: Set up nodal analysis equations.")
print("Let V_left be the voltage at the left junction and V_right be the voltage at the right junction.")
print("The top of the circuit is at 41V and the bottom is at 0V (ground).\n")

print("KCL equation at the V_left node:")
print("(V_left - 41) / Rp + (V_left - V_right) / 14 = 0")
print("\nKCL equation at the V_right node:")
print("(V_right - 41) / 11 + (V_right - V_left) / 14 + (V_right - 0) / 29 = 0")
print("-" * 30)

print("Step 3: Rearrange the equations into the standard form Ax = B and solve for the voltages.")
# Equation for V_left: V_left * (1/Rp + 1/R14) - V_right * (1/R14) = 41/Rp
# Equation for V_right: -V_left * (1/R14) + V_right * (1/R11 + 1/R14 + 1/R29) = 41/R11

# Coefficients for the matrix A
a11 = 1/Rp + 1/R14
a12 = -1/R14
a21 = -1/R14
a22 = 1/R11 + 1/R14 + 1/R29

# Vector B
b1 = V_source / Rp
b2 = V_source / R11

# Create the matrix and vector using numpy
A = np.array([[a11, a12], [a21, a22]])
B = np.array([b1, b2])

# Solve the system of linear equations
voltages = np.linalg.solve(A, B)
V_left = voltages[0]
V_right = voltages[1]

print(f"Solving the system of equations gives:")
print(f"V_left = {V_left:.4f} V")
print(f"V_right = {V_right:.4f} V")
print("-" * 30)

print("Step 4: Calculate the current I using Ohm's Law.")
print("The current I is the current flowing from V_left to V_right through the 14 Ohm resistor.")
# Calculate the current I
I = (V_left - V_right) / R14

print(f"I = (V_left - V_right) / R14")
print(f"I = ({V_left:.4f} - {V_right:.4f}) / {R14}")
print(f"I = {I:.4f} A")
print("-" * 30)

print(f"The final value of the current I is {I:.4f} A.")

final_answer = I