import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# Define the given parameters
m = 0.20  # kg, mass of the body
M = 0.80  # kg, mass of the guide
R_cm = 20   # cm, radius of the circular arcs
d_cm = 50   # cm, length of the straight section

# Convert units from cm to m
R = R_cm / 100.0
d = d_cm / 100.0

# The horizontal displacement of the guide (ΔX_M) is calculated based on the
# conservation of the system's center of mass, as there are no external horizontal forces.
# The formula derived from this principle is: ΔX_M = - (m * (R + d)) / (m + M)

# Calculate the terms for the equation
horizontal_dist_on_guide = R + d
total_mass = m + M

# Calculate the final displacement of the guide
delta_X_M = - (m * horizontal_dist_on_guide) / total_mass

# --- Outputting the solution ---
print("The formula for the horizontal displacement of the guide (ΔX_M) is:")
print("ΔX_M = - (m * (R + d)) / (m + M)\n")

print("Substituting the given values into the equation:")
# The prompt requires printing each number in the final equation
print(f"ΔX_M = - ({m} * ({R} + {d})) / ({m} + {M})")

# Show the intermediate calculation step
numerator = m * horizontal_dist_on_guide
print(f"ΔX_M = - ({numerator:.2f}) / ({total_mass:.2f})")

# Print the final result
print(f"\nThe final horizontal displacement of the guide is {delta_X_M:.2f} m.")
print("The negative sign indicates the guide moved to the left.")

# Restore stdout
sys.stdout = old_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)

# Final answer in the required format
final_answer = delta_X_M
# The problem asks for the displacement, which is a vector quantity.
# The calculated value is -0.14 m.
# We will output the numerical value.
print(f'<<<{final_answer}>>>')