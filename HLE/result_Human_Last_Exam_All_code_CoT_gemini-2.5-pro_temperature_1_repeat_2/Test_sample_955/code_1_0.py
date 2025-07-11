import cmath

# The problem analysis leads to the quadratic equation: 2*n^2 - 25*n + 50 = 0
# We need to find the integer, even solution for n.

# Coefficients of the quadratic equation an^2 + bn + c = 0
a = 2
b = -25
c = 50

# Calculate the discriminant
d = (b**2) - 4 * a * c

# Find the two solutions
sol1 = (-b - cmath.sqrt(d)) / (2 * a)
sol2 = (-b + cmath.sqrt(d)) / (2 * a)

# The problem specifies that n must be an even integer.
# We check which of the solutions meets this condition.
n = None
if sol1.imag == 0 and sol1.real == int(sol1.real) and sol1.real % 2 == 0:
    n = int(sol1.real)
elif sol2.imag == 0 and sol2.real == int(sol2.real) and sol2.real % 2 == 0:
    n = int(sol2.real)

# Now, we verify this value of n.
# Number of reachable cells R(n) = n^2/2 + 2n - 4
# Total cells = n^2
# Probability = R(n) / n^2
num_reachable = (n**2 / 2) + (2 * n) - 4
total_cells = n**2
probability = num_reachable / total_cells

print(f"The derived equation is: {a}*n^2 + ({b})*n + {c} = 0")
print(f"The solutions for n are: {sol1.real} and {sol2.real}")
print(f"The value of n that is an even integer is: {n}")
print(f"For n = {n}, the number of reachable cells is {int(num_reachable)}.")
print(f"The total number of cells is {total_cells}.")
print(f"The probability is {int(num_reachable)}/{total_cells} = {probability}, which is 66%.")
print(f"\nTherefore, the value of n is {n}.")
