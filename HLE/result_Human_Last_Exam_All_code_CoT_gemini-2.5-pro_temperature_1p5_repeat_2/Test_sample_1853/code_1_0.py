import sys

# This script calculates the gate capacitance based on Quantum Hall Effect measurements.

# --- Step 1: Define the key physical relationships ---
print("Step 1: State the fundamental equations.")
print("The carrier density 'n' can be expressed in two ways:")
print("1. From device physics: n = C_g/e * (V_bg - V_th)")
print("   where C_g is the gate capacitance we want to find.")
print("2. From the Quantum Hall effect: n = nu * (e * B / h)")
print("   where 'nu' is the Landau level filling factor.")
print("-" * 20)

# --- Step 2: Relate the change in voltage to the change in filling factor ---
print("Step 2: Find a formula for C_g by looking at the change between levels.")
print("The change in density (delta_n) for a change in voltage (delta_V_bg) is:")
print("delta_n = C_g/e * delta_V_bg")
print("This change in density corresponds to filling a new set of Landau levels:")
print("delta_n = delta_nu * (e * B / h)")
print("Equating these gives:")
print("C_g/e * delta_V_bg = delta_nu * (e * B / h)")
print("Solving for C_g:")
print("C_g = (e^2 * B / h) * (delta_nu / delta_V_bg)")
print("-" * 20)

# --- Step 3: Determine delta_V_bg and delta_nu from the problem data ---
print("Step 3: Use the values from the problem.")
# delta_V_bg
delta_V_bg_val = "3*V1 - V1"
delta_V_bg_final = "2*V1"
print(f"The gate voltages are V1, 3*V1, and 5*V1.")
print(f"The step in gate voltage is delta_V_bg = {delta_V_bg_val} = {delta_V_bg_final}.")

# delta_nu
gs = 2  # spin degeneracy
gv = 2  # valley degeneracy
delta_nu_val = gs * gv
print(f"\nThe device has spin degeneracy gs={gs} and valley degeneracy gv={gv}.")
print(f"The total degeneracy of each Landau level is g = gs * gv = {delta_nu_val}.")
print(f"The observed steps correspond to filling these full degenerate levels.")
print(f"So, the change in filling factor is delta_nu = {delta_nu_val}.")
print("-" * 20)

# --- Step 4: Substitute and calculate the final expression ---
print("Step 4: Substitute delta_nu and delta_V_bg into the formula for C_g.")
print(f"C_g = (e^2 * B / h) * ( {delta_nu_val} / ({delta_V_bg_final}) )")
print("Simplifying the expression...")

# Perform the simplification for the numeric part
numerator = 2 # from the e^2
denominator = 1 # from h*V1
numeric_factor = (numerator * delta_nu_val) / (denominator * 2) # The 2 is from 2*V1

print("\nThe final equation for the gate capacitance C_g is:")
# The prompt requests outputting each number in the final equation. The final number is 2.
final_numeric_factor = int(numeric_factor)
print(f"C_g = {final_numeric_factor} * (e^2 * B) / (h * V1)")

# For the final answer format
final_answer_expression = f"{final_numeric_factor} * e**2 * B / (h * V1)"

# Redirecting final answer to stdout
sys.stdout.write("<<<" + final_answer_expression + ">>>")