# Set the value of q
q = 997

# Calculate powers of q for convenience. Python handles large integers automatically.
q2 = q**2
q3 = q**3
q4 = q**4

# Calculate N1: the number of elements x in SU(4,q) such that x^2 = I (and x is not central)
# Formula: N1 = q^4 * (q^2 - q + 1) * (q^2 + 1)
N1 = q4 * (q2 - q + 1) * (q2 + 1)

# Calculate N2: the number of elements x in SU(4,q) such that x^2 = -I
# Formula: N2 = q^4 * (q^3 + 1)
N2 = q4 * (q3 + 1)

# The total number of involutions in PSU(4,q) is (N1 + N2) / 2
# We use integer division //
total_involutions = (N1 + N2) // 2

# Print the final result in a descriptive equation
print(f"Number of involutions = (N1 + N2) / 2 = ({N1} + {N2}) / 2 = {total_involutions}")
