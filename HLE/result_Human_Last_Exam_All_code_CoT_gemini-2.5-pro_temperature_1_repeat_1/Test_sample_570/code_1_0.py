# The minimal area is achieved by the rhombus with vertices (1,0), (0,1), (-1,0), (0,-1).
# We can calculate its area using the formula for a rhombus: Area = 1/2 * d1 * d2,
# where d1 and d2 are the lengths of the diagonals.

# The first diagonal connects (-1,0) and (1,0). Its length is:
d1 = 1 - (-1)

# The second diagonal connects (0,-1) and (0,1). Its length is:
d2 = 1 - (-1)

# Calculate the area
area = (1/2) * d1 * d2

# Print the equation and the result
print(f"The minimal area is calculated by the equation: (1/2) * {d1} * {d2} = {area}")