import math

# Define the given values from the circuit diagram
V = 41.0  # Voltage in Volts
R8 = 8.0  # Resistance in Ohms
R14 = 14.0 # Resistance in Ohms
R11 = 11.0 # Resistance in Ohms
R29 = 29.0 # Resistance in Ohms

# Step 1: Calculate the equivalent resistance of the parallel combination of R11 and R29.
R_parallel = (R11 * R29) / (R11 + R29)

# Step 2: Calculate the total resistance of the branch containing the current I.
# This branch consists of R8, R14, and the parallel part in series.
R_total_branch = R8 + R14 + R_parallel

# Step 3: Calculate the current I using Ohm's Law (I = V / R).
I = V / R_total_branch

# Output the steps of the calculation
print("Step 1: Calculate the equivalent resistance of the parallel resistors (R11 and R29).")
print(f"R_p = ({R11} * {R29}) / ({R11} + {R29})")
print(f"R_p = {R11 * R29} / {R11 + R29}")
print(f"R_p = {R_parallel:.4f} Ohms\n")

print("Step 2: Calculate the total series resistance of the branch.")
print(f"R_total = R8 + R14 + R_p")
print(f"R_total = {R8} + {R14} + {R_parallel:.4f}")
print(f"R_total = {R_total_branch:.4f} Ohms\n")

print("Step 3: Calculate the current I using Ohm's Law.")
print(f"I = V / R_total")
print(f"I = {V} / {R_total_branch:.4f}")
print(f"I = {I:.4f} A")