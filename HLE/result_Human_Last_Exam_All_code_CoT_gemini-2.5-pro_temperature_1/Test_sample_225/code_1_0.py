# The given curve is defined by the parametric equations:
# x(t) = cos(9t)
# y(t) = sin(5t)
# This is a Lissajous curve.

# For a Lissajous curve of the form (cos(at), sin(bt)), where 'a' and 'b'
# are coprime odd integers, the number of self-intersection points is given by
# the formula: (a - 1) * (b - 1) / 2.

# In this problem, we have:
a = 9
b = 5

# Check if a and b are coprime and odd.
# gcd(9, 5) = 1, so they are coprime.
# 9 and 5 are both odd.
# The conditions are satisfied, so we can apply the formula.

# Calculate the number of self-intersection points.
num_intersections = (a - 1) * (b - 1) // 2

# Print the calculation steps and the final answer.
print(f"The number of self-intersection points is calculated using the formula (a - 1) * (b - 1) / 2.")
print(f"Given a = {a} and b = {b}:")
print(f"Number of points = ({a} - 1) * ({b} - 1) / 2")
print(f"                   = ({a - 1}) * ({b - 1}) / 2")
print(f"                   = { (a - 1) * (b - 1) } / 2")
print(f"                   = {num_intersections}")