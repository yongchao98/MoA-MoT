# Define the given values from the circuit diagram
V2 = 1.0  # Voltage in Volts
R1 = 3.0  # Resistance in Ohms
R2 = 7.0  # Resistance in Ohms
R3 = 9.0  # Resistance in Ohms
R7 = 100.0 # Resistance in Ohms
R8 = 100.0 # Resistance in Ohms

# Step 1: Calculate the equivalent resistance of the parallel R1, R2, R8 combination
# R_p128 = 1 / (1/R1 + 1/R2 + 1/R8)
R_p128 = 1 / (1/R1 + 1/R2 + 1/R8)

# Step 2: Calculate the total resistance of the first major branch
# This branch has the R_p128 combination in series with R7
R_branch1 = R_p128 + R7

# The second major branch is just R3
R_branch2 = R3

# Step 3: Calculate the current in each branch
I_1 = V2 / R_branch1
I_2 = V2 / R_branch2

# Step 4: Calculate the total current by summing the currents from the parallel branches
I_total = I_1 + I_2

# Print the final equation with all the numbers
print("The total current is the sum of currents in the two main parallel branches:")
print(f"I_total = I_branch1 + I_branch2")
print(f"I_total = (V2 / R_branch1) + (V2 / R_branch2)")
print(f"R_branch1 = (1 / (1/{R1} + 1/{R2} + 1/{R8})) + {R7} = {R_branch1:.4f} Ohms")
print(f"R_branch2 = {R3} Ohms")
print(f"I_total = ({V2} V / {R_branch1:.4f} Ohm) + ({V2} V / {R_branch2:.4f} Ohm)")
print(f"I_total = {I_1:.4f} A + {I_2:.4f} A")
print(f"I_total = {I_total:.4f} A")

# <<<0.1209>>>