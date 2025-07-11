import math

# Based on an analysis using Vieta's formulas, the four roots of the
# polynomial were identified by inspecting the coefficients. The roots are:
# sqrt(14), 2*sqrt(6), sqrt(34), and 2*sqrt(11).

# To present them in the required increasing order, we can compare their squares:
# (sqrt(14))^2 = 14
# (2*sqrt(6))^2 = 4*6 = 24
# (sqrt(34))^2 = 34
# (2*sqrt(11))^2 = 4*11 = 44
# Since 14 < 24 < 34 < 44, the sorted roots are sqrt(14), 2*sqrt(6), sqrt(34), 2*sqrt(11).

# Create a list of the roots in their symbolic and numeric forms, already sorted.
sorted_roots = [
    ("sqrt(14)", math.sqrt(14)),
    ("2*sqrt(6)", 2 * math.sqrt(6)),
    ("sqrt(34)", math.sqrt(34)),
    ("2*sqrt(11)", 2 * math.sqrt(11))
]

print("The 4 roots of the equation in increasing order are:")

# The numbers that solve the equation are the roots.
# The code below prints each root's symbolic form and its numerical value.
for symbol, value in sorted_roots:
    print(f"{symbol} = {value}")
