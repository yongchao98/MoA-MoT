import math

# Set the value of q
q = 997

# The number of involutions is the sum of two terms, N_C and N_D.

# Term 1: N_C = |GU(4,q)| / |GU(2,q)|^2
# This simplifies to q**4 * (q**2 - q + 1) * (q**2 + 1)
nc_p1 = q**4
nc_p2 = q**2 - q + 1
nc_p3 = q**2 + 1
term_nc = nc_p1 * nc_p2 * nc_p3

# Term 2: N_D = |GU(4,q)| / |GL(2,q)|
# This simplifies to q**5 * (q + 1)**3 * (q**2 - q + 1) * (q**2 + 1)
nd_p1 = q**5
nd_p2 = (q + 1)**3
nd_p3 = q**2 - q + 1
nd_p4 = q**2 + 1
term_nd = nd_p1 * nd_p2 * nd_p3 * nd_p4

# Calculate the total number of involutions
total_involutions = term_nc + term_nd

# Print the final equation with all its components and the result
# The equation is N_C + N_D = total
# where N_C = nc_p1 * nc_p2 * nc_p3
# and   N_D = nd_p1 * nd_p2 * nd_p3 * nd_p4
print("The number of involutions is the sum of involutions of type C and type D.")
print("Number from Type C classes (g^2=I):")
print(f"{nc_p1} * {nc_p2} * {nc_p3} = {term_nc}")
print("\nNumber from Type D classes (g^2=-I):")
print(f"{nd_p1} * {nd_p2} * {nd_p3} * {nd_p4} = {term_nd}")
print("\nTotal number of involutions is the sum:")
print(f"{term_nc} + {term_nd} = {total_involutions}")

# For a single expression as requested:
# The format required is "each number in the final equation", which means printing the components
# of the sum.
print("\nFinal equation with all components:")
print(f"{nc_p1} * {nc_p2} * {nc_p3} + {nd_p1} * {nd_p2} * {nd_p3} * {nd_p4} = {total_involutions}")