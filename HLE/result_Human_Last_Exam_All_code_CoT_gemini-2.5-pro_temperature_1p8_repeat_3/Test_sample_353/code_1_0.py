import math

# We are calculating the angle alpha for the A(alpha) stability of BDF4.
# This angle is found to be arctan(3 * sqrt(3)).
# Here, we calculate the components that lead to this result.

# For theta = pi/3:
# R = -5/24
# I = 5 * sqrt(3) / 8
R = -5.0/24.0
I = 5.0 * math.sqrt(3) / 8.0

# The angle alpha is given by arctan(-I/R)
# where z = R + iI is on the stability boundary in the left half-plane.
tan_alpha = -I/R

# The components for the final expression arctan(3 * sqrt(3))
val1 = 3
val2_sqrt = 3
val2 = math.sqrt(val2_sqrt)

# Let's verify our computed tan_alpha matches this.
# print(f"tan(alpha) = {tan_alpha}")
# print(f"{val1} * sqrt({val2_sqrt}) = {val1*val2}")

print("The exact value of the angle alpha is given by the equation:")
print(f"alpha = arctan({val1} * sqrt({val2_sqrt}))")

# For completeness, the final value:
# final_value = math.atan(tan_alpha)
# print(f"This evaluates to approximately {math.degrees(final_value):.4f} degrees.")