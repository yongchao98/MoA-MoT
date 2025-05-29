# Given values
perimeter = 36
side1 = 7

# Solve for x
x = (perimeter - side1 - 4 - 1) / 3

# Calculate the lengths of the sides
side2 = x + 4
side3 = 2 * x + 1

# Determine the longest side
longest_side = max(side1, side2, side3)

print(longest_side)