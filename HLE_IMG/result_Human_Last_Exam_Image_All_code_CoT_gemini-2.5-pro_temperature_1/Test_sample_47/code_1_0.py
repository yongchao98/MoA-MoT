import sys

# Define the given values
V = 41.0  # Voltage in Volts
R76 = 76.0
R8 = 8.0
R14 = 14.0
R11 = 11.0
R29 = 29.0

# --- Step 1: Calculate the equivalent resistance of Branch 2 ---

# R11 and R29 are in parallel.
# Formula for parallel resistors: R_p = (R_a * R_b) / (R_a + R_b)
print("Step 1: Calculate the equivalent resistance of the parallel combination of R11 and R29.")
R_p_11_29 = (R11 * R29) / (R11 + R29)
print(f"R_p = ({R11} * {R29}) / ({R11} + {R29}) = {R_p_11_29:.4f} Ω")
print("-" * 30)

# This parallel combination is in series with R14.
# Formula for series resistors: R_s = R_a + R_b
print("Step 2: The result from Step 1 is in series with R14.")
R_s1 = R14 + R_p_11_29
print(f"R_s1 = {R14} + {R_p_11_29:.4f} = {R_s1:.4f} Ω")
print("-" * 30)

# This new combination is in series with R8.
print("Step 3: The result from Step 2 is in series with R8 to get the total resistance of Branch 2.")
R_b2 = R8 + R_s1
print(f"R_b2 = {R8} + {R_s1:.4f} = {R_b2:.4f} Ω")
print("-" * 30)

# --- Step 2: Calculate the total equivalent resistance (Req) ---

# Branch 1 (R76) is in parallel with Branch 2 (R_b2).
print("Step 4: The total equivalent resistance (Req) is R76 in parallel with Branch 2 (R_b2).")
Req = (R76 * R_b2) / (R76 + R_b2)
print(f"Req = ({R76} * {R_b2:.4f}) / ({R76} + {R_b2:.4f}) = {Req:.4f} Ω")
print("-" * 30)

# --- Step 3: Calculate the total current I ---

# Using Ohm's Law: I = V / Req
print("Step 5: Calculate the total current I using Ohm's Law (I = V / Req).")
I = V / Req
print(f"I = {V} V / {Req:.4f} Ω")
print(f"I = {I:.4f} A")
print("-" * 30)

# Final Answer
print("The final value of the current I is:")
# Use sys.stdout.write to prevent adding a newline, matching the format request.
sys.stdout.write(f"<<<{I:.4f}>>>")