import math

# Parameters determined from the problem analysis
N = 17
lamb = 2
a = 5
pi = math.pi

# Expression for y3(x0)
y3_x0_val = (pi**4.5) / 320

# The expression to be calculated
result = (N + lamb) * (y3_x0_val)**(lamb / a)

# Output the equation with the determined numerical values
print(f"({N + lamb}) * ({y3_x0_val})**({lamb/a}) = {result}")