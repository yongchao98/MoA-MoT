# Parameters from the problem statement
p = 43
n = 18
e = 3

# Derived parameters
f = n // e
q = p**f

# The number of classes is q^7 * (q-1)
# which is (p^f)^7 * (p^f - 1)
# which is p^(f*7) * (p^f - 1)
exp1 = f * 7
exp2 = f

# Calculate the final number
# Python's integers can handle arbitrary size, so direct computation is fine.
num_classes = (p**exp1) * (p**exp2 - 1)

# Output the equation with the computed numbers as requested.
print(f"The number of equivalence classes is {p}^{exp1} * ({p}^{exp2} - 1)")
print(f"This evaluates to: {num_classes}")