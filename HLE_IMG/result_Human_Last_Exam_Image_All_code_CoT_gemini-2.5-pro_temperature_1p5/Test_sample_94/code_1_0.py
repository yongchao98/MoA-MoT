import math

# Based on the reasoning above, the radius of the circles is r = 2.
r = 2

# The geometric model for the rectangle's dimensions is:
# Length X = 7 * r
# Width Y = 2 * r * (1 + sqrt(3))
X_val = 7 * r
Y_val = 2 * r * (1 + math.sqrt(3))

# Round the dimensions to the nearest integer as requested ("rounded to cm")
X_rounded = round(X_val)
Y_rounded = round(Y_val)

print("The determined radius of the circles is r = 2 cm.")
print("\nCalculating the dimensions of the outer rectangle:")
print(f"Length X = 7 * {r} = {X_val}")
print(f"Width Y = 2 * {r} * (1 + sqrt(3)) = {Y_val:.3f}")

print("\nAfter rounding to the nearest cm, the final dimensions are:")
print(f"Length X = {X_rounded} cm")
print(f"Width Y = {Y_rounded} cm")

print("\nThe final answer in the format x:y is:")
print(f"{X_rounded}:{Y_rounded}")
<<<14:11>>>