import math

# Define the coordinates
x_start = (0, 1)
y_target = (3600, 0)

# Calculate the magnitudes needed for the formula
# mag_y is the distance to the target center point
mag_y = y_target[0]

# The formula derived from the harmonic property of the probability is
# P = ln(2) / (2 * ln(3600))
log_2 = math.log(2)
log_3600 = math.log(mag_y)
probability = log_2 / (2 * log_3600)

# We print the final equation step by step, as requested.
# Note: The 'ln' denotes the natural logarithm.
print("The final probability P is calculated using the formula:")
print(f"P = ln(2) / (2 * ln({mag_y}))")
print(f"Breaking it down:")
print(f"ln(2) = {log_2}")
print(f"ln({mag_y}) = {log_3600}")
print(f"P = {log_2} / (2 * {log_3600})")
print(f"P = {log_2} / {2 * log_3600}")
print(f"P ≈ {probability}")
print("\nFinal Answer:")
# .2g format provides two significant digits
print(f"The approximate probability is {probability:.2g}")