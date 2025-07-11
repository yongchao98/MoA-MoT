import math

# From the analytical derivation, the exact solution for r_0 is (226 + 49*sqrt(2))/17.
# The numbers in this final equation are:
a = 226
b = 49
c = 2
d = 17

# As requested, here is the final equation with each number explicitly shown:
print("The final equation for the radial distance r_0 is:")
print(f"r_0 = ({a} + {b} * sqrt({c})) / {d}")

# Now, we calculate the numerical value of r_0.
r_0 = (a + b * math.sqrt(c)) / d
print("\nThe numerical value of r_0 is:")
print(r_0)