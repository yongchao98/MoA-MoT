import math

# Step 1 & 2: Define relationships based on geometry and the 4% area clue.
# The relationship 11*gw = 9*r was derived from:
# y = 6*r
# x = 6*r + gw
# gh = 2*r
# (gw * gh) / (x * y) = 1/25 (which is 4%)
# (gw * 2*r) / ((6*r + gw) * 6*r) = 1/25
# gw / (3 * (6*r + gw)) = 1/25
# 25*gw = 18*r + 3*gw
# 22*gw = 18*r
# 11*gw = 9*r

# Step 3: Find the smallest integer solutions for r and gw.
# Since r and gw are integers, the simplest solution to 11*gw = 9*r is:
gw = 9
r = 11
print(f"From the relationship 11*gw = 9*r, and knowing r and gw are integers:")
print(f"The smallest integer radius of a circle (r) is {r} cm.")
print(f"The smallest integer width of a green rectangle (gw) is {gw} cm.\n")

# Step 4: Calculate the final dimensions of the outer rectangle.
x = 6 * r + gw
y = 6 * r

print(f"Calculating the dimensions of the outer rectangle:")
print(f"Length x = 6 * r + gw = 6 * {r} + {gw} = {x} cm.")
print(f"Width y = 6 * r = 6 * {r} = {y} cm.\n")

# Step 5: State the final answer in the requested format.
print("The final answer for the size of the outer rectangle (x:y) is:")
print(f"{x}:{y}")

# Verification step (optional, to check against Monte Carlo data)
# Area_total = x * y
# Area_circles = 9 * math.pi * r**2
# ratio = Area_circles / Area_total
# print(f"\nVerification: The model's circle area ratio is {ratio:.3f}, compared to the Monte Carlo average of ~0.735.")
