# Define the values of the components from the circuit diagram
V2 = 1.0  # Voltage in Volts
R1 = 3.0  # Resistance in Ohms
R2 = 7.0  # Resistance in Ohms
R3 = 9.0  # Resistance in Ohms
R7 = 100.0 # Resistance in Ohms
R8 = 100.0 # Resistance in Ohms

# Step 1: Calculate the equivalent resistance of the parallel group of R1, R2, and R8.
print("Step 1: Resistors R1, R2, and R8 are in parallel.")
print(f"1/Rp1 = 1/R1 + 1/R2 + 1/R8 = 1/{R1} + 1/{R2} + 1/{R8}")
inv_Rp1 = 1/R1 + 1/R2 + 1/R8
Rp1 = 1 / inv_Rp1
print(f"Rp1 = {Rp1:.4f} Ω\n")

# Step 2: This parallel group (Rp1) is in series with R7.
print("Step 2: The parallel group (Rp1) is in series with R7.")
print(f"R_branch2 = Rp1 + R7 = {Rp1:.4f} + {R7}")
R_branch2 = Rp1 + R7
print(f"R_branch2 = {R_branch2:.4f} Ω\n")

# Step 3: The branch containing R_branch2 is in parallel with R3. Calculate the total equivalent resistance (Req).
print("Step 3: This entire branch (R_branch2) is in parallel with R3.")
print(f"1/Req = 1/R_branch2 + 1/R3 = 1/{R_branch2:.4f} + 1/{R3}")
inv_Req = 1/R_branch2 + 1/R3
Req = 1 / inv_Req
print(f"Req = {Req:.4f} Ω\n")

# Step 4: Use Ohm's Law to find the total current from the voltage source V2.
print("Step 4: Use Ohm's Law to find the total current (I_total).")
print(f"I_total = V2 / Req = {V2} V / {Req:.4f} Ω")
I_total = V2 / Req
print(f"The total current flowing through the circuit is: {I_total:.4f} A")