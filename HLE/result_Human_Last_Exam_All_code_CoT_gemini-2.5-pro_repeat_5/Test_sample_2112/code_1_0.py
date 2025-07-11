import math

# The value r_0 is derived from the analysis of the functional equation.
# The exact solution for r_0 is (226 + 49 * sqrt(2)) / 17.
a = 226
b = 49
c = 2
d = 17

# We calculate the numerical value of r_0.
r_0 = (a + b * math.sqrt(c)) / d

print("The radial distance r_0 is found from the equation:")
print(f"r_0 = ({a} + {b} * sqrt({c})) / {d}")
print(f"The calculated value of the radial distance r_0 is: {r_0}")