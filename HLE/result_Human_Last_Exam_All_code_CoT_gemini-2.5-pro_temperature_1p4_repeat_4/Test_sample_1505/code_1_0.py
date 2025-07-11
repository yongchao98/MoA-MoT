# A positive integer n for the calculation.
# You can change this value to any positive integer.
n = 10

print(f"Calculating the sum for n = {n}\n")

# The approximation formula derived from the Euler-Maclaurin expansion is:
# S_n ≈ (n^2 / 2) + (1 / 120) + (1 / (252 * n))
# The numbers in this equation are 2, 1, 120, 1, 252.

# Coefficients and constants from the formula
c1_den = 2
c2_num = 1
c2_den = 120
c3_num = 1
c3_den = 252

# Calculate the approximation for the given n
approximation = (n**2 / c1_den) + (c2_num / c2_den) + (c3_num / (c3_den * n))

# Output the final equation with each number explicitly shown, as requested.
print("The final equation for the approximation is:")
print(f"({n}**2 / {c1_den}) + ({c2_num} / {c2_den}) + ({c3_num} / ({c3_den} * {n}))")

# Print the final numerical result
print(f"\nThe calculated value is: {approximation}")