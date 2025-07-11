#
# Plan:
# 1. Define the physical parameters for the transistor with example values.
# 2. Use the formula for displacement field (D) based on the charge-control model.
#    The total displacement field is the sum of the fields from the top and back gates.
#    Formula: D = C_tg * V_tg + C_bg * V_bg
# 3. Calculate the result.
# 4. Print the steps of the calculation, including the formula with numerical values.
#

# Step 1: Define the given parameters with example values
# C_tg: Top gate capacitance per unit area (in F/m^2)
C_tg = 3.5e-3
# V_tg: Top gate voltage (in V)
V_tg = 1.0
# C_bg: Back gate capacitance per unit area (in F/m^2)
C_bg = 0.5e-3
# V_bg: Back gate voltage (in V)
V_bg = 0.5

# The dielectric constant of the transistor, epsilon_s, is not needed for this calculation
# as we are finding the net displacement field that induces charge in the grounded channel.

# Step 2: Calculate the total displacement field (D)
D = C_tg * V_tg + C_bg * V_bg

# Step 3 & 4: Print the equation and the final result
print("The formula for the total displacement field (D) that induces charge in the channel is:")
print("D = C_tg * V_tg + C_bg * V_bg")
print("\nSubstituting the given values into the equation:")
# The final code outputs each number in the final equation as requested.
print(f"D = {C_tg} F/m^2 * {V_tg} V + {C_bg} F/m^2 * {V_bg} V")

# The final result
print(f"\nThe calculated total displacement field is:")
print(f"D = {D} C/m^2")