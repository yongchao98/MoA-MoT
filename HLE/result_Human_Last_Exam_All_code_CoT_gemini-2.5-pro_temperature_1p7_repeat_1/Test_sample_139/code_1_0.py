# Define the resistor values found through analysis
R1 = 10  # ohms
R2 = 6   # ohms
R3 = 13  # ohms, which is a prime number

# Define the voltage across R3 when R2 fails
V_fail = 26  # volts

# --- Calculation for the Parallel Circuit Model ---
# This model yields the maximum possible current.

# Step 1: When R2 fails, R1 and R3 are in parallel with the current source.
# Calculate the source current based on this failure condition.
# I_src = V_fail / R_eq_fail = V_fail / (1/R1 + 1/R3)^-1 = V_fail * (1/R1 + 1/R3)
I_src = V_fail * (1/R1 + 1/R3)

# Step 2: When R2 is intact, all three resistors are in parallel.
# Calculate the equivalent resistance of the intact circuit.
R_eq_intact = 1 / (1/R1 + 1/R2 + 1/R3)

# Step 3: Calculate the voltage across the parallel bank in the intact circuit.
V_intact = I_src * R_eq_intact

# Step 4: Calculate the current through R3 in the intact circuit.
I3_intact = V_intact / R3

# Output the results and the equation for the final answer
# The user wants to see the numbers in the final equation.
# I3_intact = (V_fail * (1/R1 + 1/R3)) * (1 / (1/R1 + 1/R2 + 1/R3)) / R3
# Let's show the final calculation with numbers
numerator_I_src = V_fail * (R1 + R3)
denominator_I_src = R1 * R3

# R_eq_intact calculation:
# Numerator is R1*R2*R3
# Denominator is R1*R2 + R1*R3 + R2*R3
num_Req = R1 * R2 * R3
den_Req = R1 * R2 + R1 * R3 + R2 * R3

# V_intact = (numerator_I_src / denominator_I_src) * (num_Req / den_Req)
# I3_intact = V_intact / R3
final_numerator = (V_fail * (R1 + R3) / (R1 * R3)) * (R1 * R2 * R3 / (R1*R2 + R1*R3 + R2*R3)) / R3
final_numerator_simplified = V_fail * (R1 + R3) * R2

final_denominator_simplified = (R1 * R2 + R1 * R3 + R2 * R3)

print(f"The resistor values that maximize the current are R1 = {R1} ohms, R2 = {R2} ohms, and R3 = {R3} ohms.")
print("The circuit configuration yielding the maximum current is the parallel model.")
print("The calculation for the current through R3 is:")
print(f"I3 = (V_fail * (R1 + R3) * R2) / (R1*R2 + R1*R3 + R2*R3)")
print(f"I3 = ({V_fail} * ({R1} + {R3}) * {R2}) / ({R1}*{R2} + {R1}*{R3} + {R2}*{R3})")
print(f"I3 = ({V_fail} * {R1+R3} * {R2}) / ({R1*R2} + {R1*R3} + {R2*R3})")
print(f"I3 = ({final_numerator_simplified}) / ({final_denominator_simplified})")
print(f"The maximum possible current through R3 is {I3_intact} Amperes.")