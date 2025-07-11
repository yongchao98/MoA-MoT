# The prime q for PSU(4, q) is given
q = 997

# Calculate N1, the size of the first conjugacy class of involutions
# N1 = q^4 * (q^2 + 1) * (q^2 - q + 1)
n1 = q**4 * (q**2 + 1) * (q**2 - q + 1)

# Calculate N2, the size of the second conjugacy class of involutions
# N2 = q^5 * (q + 1)^2 * (q^2 + 1) * (q^2 - q + 1)
n2 = q**5 * (q + 1)**2 * (q**2 + 1) * (q**2 - q + 1)

# The total number of involutions is the sum of the sizes of the two classes
total_involutions = n1 + n2

# Print the final equation as requested
print(f"{n1} + {n2} = {total_involutions}")

# Return the final numerical answer in the specified format
# print(f"<<<{total_involutions}>>>")