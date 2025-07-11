# Set the parameters based on the plan
# L: Number of loops
# n: Number of lines meeting at a vertex (for phi-4 theory)
# E: Number of external lines (for a vacuum bubble diagram)
L = 2
n = 4
E = 0

# The two fundamental equations for a Feynman diagram are:
# 1) L = I - V + 1
# 2) n * V = 2 * I + E
#
# From (1), we can express the number of internal lines (I) as:
# I = L + V - 1
#
# Substitute this expression for I into equation (2):
# n * V = 2 * (L + V - 1) + E
#
# Now, we solve for V (the number of vertices):
# n * V = 2*L + 2*V - 2 + E
# n*V - 2*V = 2*L - 2 + E
# V * (n - 2) = 2*L - 2 + E
# V = (2*L - 2 + E) / (n - 2)

# Now, we calculate the numerator and denominator separately
# to show the steps of the calculation.
numerator = 2 * L - 2 + E
denominator = n - 2
min_vertices = numerator / denominator

# Print the final calculation step-by-step
print("To find the minimum number of vertices (V), we use the formula:")
print("V = (2*L - 2 + E) / (n - 2)\n")
print(f"Given L = {L} (two loops), n = {n} (phi-4 theory), and E = {E} (vacuum diagram):")
print(f"V = (2*{L} - 2 + {E}) / ({n} - 2)")
print(f"V = ({2*L - 2 + E}) / ({n-2})")
print(f"V = {int(min_vertices)}")
print(f"\nThe minimum number of vertices in a two-loop Feynman diagram is {int(min_vertices)}.")
