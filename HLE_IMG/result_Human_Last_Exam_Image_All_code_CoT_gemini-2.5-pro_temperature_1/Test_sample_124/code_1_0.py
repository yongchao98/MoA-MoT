import sys

# Define the component values from the circuit diagram
V2 = 1.0  # Voltage in Volts
R1 = 3.0  # Resistance in Ohms
R2 = 7.0  # Resistance in Ohms
R3 = 9.0  # Resistance in Ohms
R7 = 100.0 # Resistance in Ohms
R8 = 100.0 # Resistance in Ohms

# Step 1: Explain the circuit configuration.
print("The circuit consists of four parallel branches connected to the 1V source:")
print(f"1. Resistor R1 = {R1} Ω")
print(f"2. Resistor R2 = {R2} Ω")
print(f"3. Resistor R3 = {R3} Ω")
print(f"4. A series combination of R7 = {R7} Ω and R8 = {R8} Ω")
print("-" * 30)

# Step 2: Calculate the resistance of the series branch.
R_series_78 = R7 + R8
print(f"The resistance of the series branch is R7 + R8 = {R7} + {R8} = {R_series_78} Ω")
print("-" * 30)

# Step 3: Calculate the current through each parallel branch.
I1 = V2 / R1
I2 = V2 / R2
I3 = V2 / R3
I4 = V2 / R_series_78

# Step 4: Calculate the total current by summing the currents of the parallel branches.
total_current = I1 + I2 + I3 + I4

# Step 5: Print the final equation and the result.
print("The total current is the sum of the currents in each branch:")
# Using sys.stdout.write to avoid adding newlines automatically, for equation formatting
sys.stdout.write(f"I_total = (V2 / R1) + (V2 / R2) + (V2 / R3) + (V2 / (R7 + R8))\n")
sys.stdout.write(f"I_total = ({V2} / {R1}) + ({V2} / {R2}) + ({V2} / {R3}) + ({V2} / ({R7} + {R8}))\n")
sys.stdout.write(f"I_total = {I1:.4f} A + {I2:.4f} A + {I3:.4f} A + {I4:.4f} A\n")
print("-" * 30)
print(f"The total current flowing through the circuit is: {total_current:.4f} A")