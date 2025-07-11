# The problem reduces to calculating the number of conjugacy classes
# of the extraspecial group S = 3^(1+2)_+.

# Parameters for the extraspecial group p^(1+2n)_+
p = 3
n = 1

# The number of conjugacy classes of an extraspecial p-group of order p^(2n+1)
# is given by the formula: p + (p^(2n) - 1).
num_classes = p + (p**(2 * n) - 1)

# We print the final calculation step-by-step as requested.
# First part of the equation: p
part1 = p
# Second part of the equation: p^(2n)
part2 = p**(2 * n)
# Third part of the equation: - 1
part3 = -1

# The expression inside the parenthesis
inside_parenthesis = part2 + part3

print(f"The number of blocks is found by calculating the number of conjugacy classes of S.")
print(f"The formula is p + (p^(2n) - 1). With p={p} and n={n}, the calculation is:")
print(f"{part1} + ({p}^(2*{n}) - 1) = {part1} + ({part2} - 1) = {part1} + {inside_parenthesis} = {num_classes}")
print(f"\nThus, kG has {num_classes} blocks.")
