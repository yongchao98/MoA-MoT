import sys
import io

# Backup the standard output
stdout_backup = sys.stdout
# Redirect standard output to a string
sys.stdout = captured_output = io.StringIO()

# Parameters from the problem
m = 0.20  # kg, mass of the body
M = 0.80  # kg, mass of the guide
R_cm = 20  # cm, radius of the circular arcs
d_cm = 50  # cm, length of the straight section

# Convert units from cm to m
R = R_cm / 100.0
d = d_cm / 100.0

# The horizontal displacement of the mass relative to the guide (Δx_m_rel)
# is the sum of the horizontal projections of the path segments.
# Δx_m_rel = R (first arc) + d (straight section) + R (second arc)
delta_x_m_rel = 2 * R + d

# The formula for the guide's displacement (Δx_g) derived from the
# conservation of the center of mass position is:
# Δx_g = - (m / (m + M)) * Δx_m_rel
# We will now substitute the values into this formula.

print("The formula for the guide's displacement (Δx_g) is:")
print("Δx_g = - (m / (m + M)) * (2*R + d)")
print("\nPlugging in the values:")
print(f"Δx_g = - ({m} kg / ({m} kg + {M} kg)) * (2 * {R} m + {d} m)")

# Calculate the intermediate steps to show the process clearly
m_plus_M = m + M
two_R_plus_d = 2 * R + d
print(f"Δx_g = - ({m:.2f} / {m_plus_M:.2f}) * ({two_R_plus_d:.2f})")

ratio = m / m_plus_M
print(f"Δx_g = - {ratio:.2f} * {two_R_plus_d:.2f}")

# Calculate the final result
delta_x_g = - (m / (m + M)) * (2 * R + d)

# Print the final result
print(f"The final calculated displacement of the guide is:")
print(f"Δx_g = {delta_x_g:.2f} m")
print("\nThe negative sign indicates the guide moves to the left, opposite to the mass's relative displacement.")

# Restore the standard output and print the captured string
sys.stdout = stdout_backup
output_str = captured_output.getvalue()
print(output_str)

print("<<<" + f"{delta_x_g:.2f}" + ">>>")