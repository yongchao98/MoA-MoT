import math

# The four roots are identified from the polynomial's coefficients as
# sqrt(14), 2*sqrt(6), sqrt(34), and 2*sqrt(11).
# The final equation in factored form is (X - r_1)(X - r_2)(X - r_3)(X - r_4) = 0.
# The numbers in this equation are the roots r_1, r_2, r_3, r_4.
# This code calculates these four numbers and prints them in increasing order.

roots = [
    math.sqrt(14),
    2 * math.sqrt(6),
    math.sqrt(34),
    2 * math.sqrt(11)
]

# Sort the roots into increasing order
roots.sort()

# Output each number (root) in the final equation.
# The roots are printed in increasing order as requested.
print("The 4 roots of the equation in increasing order are:")
print(roots[0])
print(roots[1])
print(roots[2])
print(roots[3])