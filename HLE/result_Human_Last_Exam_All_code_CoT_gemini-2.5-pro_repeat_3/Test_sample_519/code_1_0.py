# This script formats the properties of the three given CFGs and prints the result.

# Properties for X1: Hilb^11(A^3)
# The dimension is calculated as n * d = 3 * 11 = 33.
profile_x1 = "[S, s, 33]"

# Properties for X2: [(A^4 \ V(xy-zw)) / C*]
# The dimension is calculated as dim(space) - dim(group) = 4 - 1 = 3.
profile_x2 = "[DM, s, irr, 3]"

# Properties for X3: The Picard stack of a genus 7 curve.
# The dimension is the genus of the curve, which is 7.
profile_x3 = "[A, s, 7]"

# Combine the profiles into a single string, separated by spaces.
final_output = f"{profile_x1} {profile_x2} {profile_x3}"

print(final_output)