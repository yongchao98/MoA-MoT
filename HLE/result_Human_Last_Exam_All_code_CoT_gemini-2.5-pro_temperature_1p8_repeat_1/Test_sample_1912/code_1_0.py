# The input numbers for the final equation
a = 6
b = 12

# The discovered pattern is a * b + a - (b / a).
# In Python, we use // for integer division.
result = a * b + a - (b // a)

# Print out the components of the final equation as requested
# The final equation is 6 + 12 = result
# The pattern used is a * b + a - (b / a)
# Let's show the full calculation.
print(f"The pattern is: result = a * b + a - (b / a)")
print(f"For the final equation, 6 + 12 = ?:")
print(f"a = {a}, b = {b}")
print(f"Result = {a} * {b} + {a} - ({b} / {a})")
print(f"Result = {a*b} + {a} - {b//a}")
print(f"Result = {a*b + a} - {b//a}")
print(f"Result = {result}")
print(f"So, the final equation is:")
print(f"6 + 12 = {result}")