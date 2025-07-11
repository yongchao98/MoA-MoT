import math

# Define the constants from the derived formula for y(0)
# y(0) = (a/b) * (c/d)^(e/f)
a = 3
b = 4
c = 9
d = 4
e = 3
f = 5

# Calculate the value of y(0)
# The formula is y(0) = (3/4) * (9/4)^(3/5)
result = (a / b) * ((c / d)**(e / f))

# Print the final equation with its components
print(f"The final equation for y(0) is: y(0) = ({a}/{b}) * ({c}/{d})^({e}/{f})")

# Print the calculated value of y(0)
print(f"The membrane's deflection at x = 0 is y(0) = {result}")
