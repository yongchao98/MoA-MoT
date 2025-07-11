import math

# Coefficients of the quadratic equation r^2 - 97r - 61 = 0
a = 1
b = -97
c = -61

# Calculate the discriminant
delta = b**2 - 4*a*c

# Find the two roots
# We need the positive root which corresponds to the '+' sign
r0 = (-b + math.sqrt(delta)) / (2*a)

print(f"The radial distance r_0 is determined by the quadratic equation:")
print(f"{a} * r_0^2 + ({b}) * r_0 + ({c}) = 0")
print("\nThe full equation for r_0 is:")
print(f"r_0 = ( -({b}) + sqrt( ({b})^2 - 4*({a})*({c}) ) ) / ( 2*({a}) )")
print("\nSolving for r_0 > 15, we get:")
print(f"r_0 = {r0}")
