import sys
# Redirect stdout to a file to prevent unwanted output
original_stdout = sys.stdout
sys.stdout = open('output.txt', 'w')
import numpy as np

# Step 1: Define component values from the circuit diagram.
V = 1.0   # Voltage in Volts
R1 = 3.0  # Resistance in Ohms
R2 = 7.0  # Resistance in Ohms
R3 = 9.0  # Resistance in Ohms
R7 = 100.0# Resistance in Ohms
R8 = 100.0# Resistance in Ohms

# Step 2: The circuit analysis shows that R1, R2, R3, and the series
# combination of (R7+R8) are all in parallel with the voltage source.
# The total current is the sum of the currents through each parallel branch.

# First, calculate the combined resistance of the R7-R8 series branch.
R_7_8 = R7 + R8

# Step 3: Calculate the current through each parallel branch using Ohm's Law (I = V/R).
I1 = V / R1
I2 = V / R2
I3 = V / R3
I_7_8 = V / R_7_8

# Step 4: Calculate the total current by summing the individual branch currents.
I_total = I1 + I2 + I3 + I_7_8

# Step 5: Print the final calculation, showing the equation with all the numbers.
# We are calculating: I_total = V/R1 + V/R2 + V/R3 + V/(R7+R8)
sys.stdout = original_stdout
print("The total current (I_total) is found by summing the currents in each parallel branch.")
print(f"I_total = V/R1 + V/R2 + V/R3 + V/(R7+R8)")
print(f"I_total = {V}/{R1} + {V}/{R2} + {V}/{R3} + {V}/({R7}+{R8})")
print(f"I_total = {I_total:.4f} A")
print(f"<<<{I_total:.4f}>>>")