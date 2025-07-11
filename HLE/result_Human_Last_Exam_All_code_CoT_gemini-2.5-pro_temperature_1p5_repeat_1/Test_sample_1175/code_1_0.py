import math

# The final expression is of the form:
# (a/b) * base^(p1/p2) + (c/d)

a = 3
b = 2
base = 10
p1 = 10
p2 = 3
c = 37
d = 4

# Print the final symbolic expression
print(f"The final expression is: ({a}/{b}) * {base}^({p1}/{p2}) + ({c}/{d})")

# Calculate the numerical value
result = (a / b) * math.pow(base, p1 / p2) + (c / d)

# Print the final numerical result
print(f"The numerical value is: {result}")