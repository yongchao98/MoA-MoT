import math

# Step 1: Define the problem based on geometric interpretation.
# The problem asks for the minimal area of a compact set C such that for any direction,
# the length of the longest chord in that direction is at least w = 1/2.

# Step 2: State the relevant theorem.
# According to a theorem by Fujiwara, the convex set with the smallest area that satisfies
# this property is an equilateral triangle of height w.

# Step 3: Define the given parameter.
# The minimum required length for the longest chord is:
h = 0.5

# Step 4: Define the formula for the area of an equilateral triangle from its height.
# The area 'A' of an equilateral triangle with height 'h' is A = h^2 / sqrt(3).
print("The minimal area is achieved by an equilateral triangle of height h = 1/2.")
print("The formula for the area (A) of such a triangle is: A = h² / √3")
print("")

# Step 5: Calculate the area and display the equation with numbers.
sqrt3_val = math.sqrt(3)
area = h**2 / sqrt3_val

print("Plugging in the value h = 0.5:")
print(f"A = ({h})² / √3")
print(f"A = {h**2} / {sqrt3_val}")
print(f"A = {area}")
